/**
 * @file bucket.cpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <tascel.h>

#include "constants.h"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"

#define USE_ARMCI 1
#if USE_ARMCI
/* ANL's armci does not extern "C" inside the header */
extern "C" {
#include <armci.h>
}
#endif

using tascel::LockGuard;
using tascel::PthreadMutex;

namespace pgraph {

/** power function for size_t */
static size_t powz(size_t base, size_t n)
{
    size_t p = 1;

    for(/*empty*/; n > 0; --n) {
        assert(p < SIZE_MAX/base);
        p *= base;
    }

    return p;
}


/**
 * This function can be optimized if higher performance is required.
 *
 * @param[in] kmer address of of k-mer string
 * @param[in] k slide window size
 * @return TODO todo
 */
static inline size_t entry_index(const char *kmer, int k)
{
    int i;
    size_t value = 0;

    for (i = 0; i < k; ++i) {
        value = value * SIGMA + (kmer[i] - 'A');
    }

    return value;
}


static bool SuffixBucketIndexCompare(const Suffix &i, const Suffix &j)
{
    return i.bid < j.bid;
}


ostream& operator<<(ostream &os, const Suffix &suffix)
{
    return os << "(" << suffix.sid
        << "," << suffix.pid
        << "," << suffix.bid
        << ")";
}


SuffixBuckets::SuffixBuckets(SequenceDatabase *sequences,
                             const Parameters &param,
                             MPI_Comm comm,
                             bool dumb_split)
    :   sequences(sequences)
    ,   param(param)
    ,   suffixes(NULL)
    ,   suffixes_size(0)
    ,   buckets(NULL)
    ,   buckets_size(0)
    ,   first_bucket(0)
    ,   last_bucket(0)
    ,   dumb_split(dumb_split)
{
    size_t n_suffixes = 0;
    size_t n_buckets = 0;
    size_t suffix_index = 0;
    size_t n_seq = 0;
    size_t remainder = 0;
    size_t start = 0;
    size_t stop = 0;
    int ierr = 0;

    ierr = MPI_Comm_rank(comm, &comm_rank);
    assert(MPI_SUCCESS == ierr);
    ierr = MPI_Comm_size(comm, &comm_size);
    assert(MPI_SUCCESS == ierr);

    /* allocate buckets */
    n_buckets = powz(SIGMA, param.window_size);
    buckets = new Bucket[n_buckets];
    bucket_size.assign(n_buckets, size_t(0));
    bucket_owner.assign(n_buckets, size_t(0));
    bucket_offset.assign(n_buckets, size_t(0));

    /* allocate suffixes */
    n_suffixes = sequences->get_global_size()
        - sequences->get_global_count() * param.window_size;

    mpix_print_zero("n_buckets", n_buckets, comm);
    mpix_print_zero("n_suffixes", n_suffixes, comm);

    /* each MPI rank gets a contiguous range of sequences to bucket */
    n_seq = sequences->get_global_count() / comm_size;
    remainder = sequences->get_global_count() % comm_size;
    start = n_seq * comm_rank;
    stop = n_seq * (comm_rank+1);
    if (comm_rank < remainder) {
        start += comm_rank;
        stop += comm_rank+1;
    }
    else {
        start += remainder;
        stop += remainder;
    }
    if (stop > sequences->get_global_count()) {
        stop = sequences->get_global_count();
    }
#if DEBUG
    mpix_print_sync("n_seq", n_seq, comm);
    mpix_print_sync("remainder", remainder, comm);
    mpix_print_sync("start", start, comm);
    mpix_print_sync("stop", stop, comm);
#endif

    size_t initial_suffixes_size = 0;
    for (size_t i = start; i < stop; ++i) {
        Sequence &sequence = sequences->get_sequence(i);
        initial_suffixes_size += sequence.get_sequence_length();
    }
    vector<Suffix> initial_suffixes(initial_suffixes_size);

    /* slide k-mers for every sequence and bucket them */
    for (size_t i = start; i < stop; ++i) {
        Sequence &sequence = sequences->get_sequence(i);
        size_t stop_index = 0;
        const char *sequence_data = NULL;
        size_t sequence_length = 0;

        sequence.get_sequence(sequence_data, sequence_length);
        stop_index = sequence_length - param.window_size - 1;

        if (sequence_length <= param.window_size) continue;

        for (size_t j = 0; j <= stop_index; ++j) {
            size_t bucket_index = entry_index(
                    sequence_data + j, param.window_size);
            assert(bucket_index < n_buckets);
            /* prefixed in the suffix list for the given bucket */
            initial_suffixes[suffix_index].sid = i;
            initial_suffixes[suffix_index].pid = j;
            initial_suffixes[suffix_index].bid = bucket_index;
            initial_suffixes[suffix_index].next = NULL;
            bucket_size[bucket_index]++;
            buckets[bucket_index].size++;
            suffix_index++;
        }
    }
    initial_suffixes.resize(suffix_index);

    mpix_allreduce(bucket_size, MPI_SUM, comm);
#if DEBUG
    {
        size_t count=0;
        for (size_t i=0; i<n_buckets; ++i) {
            count += bucket_size[i];
        }
        mpix_print_sync("count", count, comm);
    }
#endif
    //assert(count == n_suffixes); may not be true if a seq is < window size

#if DEBUG
    mpix_print_sync("suffix_index", suffix_index, comm);
#endif

    n_seq = sequences->get_global_count() / comm_size;
    remainder = sequences->get_global_count() % comm_size;
    start = n_seq * comm_rank;
    stop = n_seq * (comm_rank+1);
    if (comm_rank < remainder) {
        start += comm_rank;
        stop += comm_rank+1;
    }
    else {
        start += remainder;
        stop += remainder;
    }
    if (dumb_split) {
        size_t even_split = n_buckets / comm_size;
        size_t remainder = n_buckets % comm_size;
        size_t start = even_split * comm_rank;
        size_t stop = even_split * (comm_rank+1);
        if (comm_rank < remainder) {
            start += comm_rank;
            stop += comm_rank+1;
        }
        else {
            start += remainder;
            stop += remainder;
        }
        mpix_print_sync("n_buckets", n_buckets, comm);
        mpix_print_sync("comm_size", comm_size, comm);
        mpix_print_sync("even_split", even_split, comm);
        mpix_print_sync("remainder", remainder, comm);
        mpix_print_sync("start", start, comm);
        mpix_print_sync("stop", stop, comm);
        for (size_t i=start; i<stop; ++i) {
            bucket_owner[i] = comm_rank;
        }
        mpix_allreduce(bucket_owner, MPI_SUM, comm);
    }
    else {
        size_t even_split = n_suffixes / comm_size;
        size_t rank = 0;
        size_t count = 0;
        for (size_t i=0; i<n_buckets; ++i) {
            count += bucket_size[i];
            if (count > even_split) {
                count = 0;
                ++rank;
                if (rank >= comm_size) {
                    rank = comm_size-1;
                }
                if (rank == comm_rank) {
                    first_bucket = i;
                }
            }
            bucket_owner[i] = rank;
        }
    }

    vector<int> amount_to_send(comm_size, 0);
    vector<int> amount_to_recv(comm_size, 0);
    for (size_t i=0; i<n_buckets; ++i) {
        amount_to_send[bucket_owner[i]] += int(buckets[i].size);
        buckets[i].size = 0; /* reset for later */
    }
#if DEBUG
    mpix_print_sync("amount_to_send", vec_to_string(amount_to_send), comm);
#endif
    mpix_alltoall(amount_to_send, amount_to_recv, comm);
#if DEBUG
    mpix_print_sync("amount_to_recv", vec_to_string(amount_to_recv), comm);
#endif

    int total_amount_to_recv = 0;
    for (int i=0; i<comm_size; ++i) {
        total_amount_to_recv += amount_to_recv[i];
    }

    /* We are preparing for the all to all, so we sort the suffixes. Yes, this
     * invalidates any of their 'next' pointers but these will be rebuilt after
     * the all to all. */
    std::sort(initial_suffixes.begin(),
              initial_suffixes.end(),
              SuffixBucketIndexCompare);
#if USE_ARMCI
    bucket_address.assign(comm_size, NULL);
    (void)ARMCI_Malloc((void**)&bucket_address[0],
            total_amount_to_recv*sizeof(Suffix));
    suffixes = bucket_address[comm_rank];
#else
    suffixes = new Suffix[total_amount_to_recv];
#endif
    vector<int> send_displacements(comm_size, 0);
    vector<int> recv_displacements(comm_size, 0);
    for (size_t i=1; i<comm_size; ++i) {
        send_displacements[i] = send_displacements[i-1] + amount_to_send[i-1];
        recv_displacements[i] = recv_displacements[i-1] + amount_to_recv[i-1];
    }
#if DEBUG
    mpix_print_sync("send_displacements", vec_to_string(send_displacements), comm);
    mpix_print_sync("recv_displacements", vec_to_string(recv_displacements), comm);
#endif

    /* need to alltoallv the buckets to the owning processes
     * probably need to add a "bucket_id" field to the Suffix class so that
     * after the data is sent it can be properly recovered to the correct
     * bucket linked lists */
    MPI_Datatype SuffixType;
    MPI_Datatype type[4] = {
        mpix_get_mpi_datatype(initial_suffixes[0].sid),
        mpix_get_mpi_datatype(initial_suffixes[0].pid),
        mpix_get_mpi_datatype(initial_suffixes[0].bid),
        MPI_AINT
    };
    int blocklen[4] = {1,1,1,1};
    MPI_Aint disp[4];
    disp[0] = MPI_Aint(&initial_suffixes[0].sid) - MPI_Aint(&initial_suffixes[0]);
    disp[1] = MPI_Aint(&initial_suffixes[0].pid) - MPI_Aint(&initial_suffixes[0]);
    disp[2] = MPI_Aint(&initial_suffixes[0].bid) - MPI_Aint(&initial_suffixes[0]);
    disp[3] = MPI_Aint(&initial_suffixes[0].next)- MPI_Aint(&initial_suffixes[0]);
    ierr = MPI_Type_create_struct(4, blocklen, disp, type, &SuffixType);
    assert(MPI_SUCCESS == ierr);
    ierr = MPI_Type_commit(&SuffixType);
    assert(MPI_SUCCESS == ierr);
    //mpix_print_sync("initial_suffixes", vec_to_string(initial_suffixes), comm);
    ierr = MPI_Alltoallv(
            &initial_suffixes[0], &amount_to_send[0],
            &send_displacements[0], SuffixType,
            &suffixes[0], &amount_to_recv[0],
            &recv_displacements[0], SuffixType, comm);
    assert(MPI_SUCCESS == ierr);
    //mpix_print_sync("suffixes", arr_to_string(suffixes, total_amount_to_recv), comm);

#if USE_ARMCI
    std::sort(suffixes,
              suffixes+total_amount_to_recv,
              SuffixBucketIndexCompare);
    /* The suffixes contains sorted suffixes based on the buckets they belong
     * to. That means we can simply update the 'next' links! */
    size_t last_id = sequences->get_global_count();
    bucket_size_total = 0;
    for (size_t i=0; i<total_amount_to_recv-1; ++i) {
        assert(suffixes[i].sid < last_id);
        assert(suffixes[i].bid < n_buckets);
        assert(suffixes[i].bid <= suffixes[i+1].bid);
        assert(suffixes[i].next == NULL);
        if (suffixes[i].bid == suffixes[i+1].bid) {
            suffixes[i].next = &suffixes[i+1];
        }
        else {
            suffixes[i].next = NULL;
        }
        if (buckets[suffixes[i].bid].suffixes == NULL) {
            buckets[suffixes[i].bid].suffixes = &suffixes[i];
            bucket_offset[suffixes[i].bid] = i;
        }
        buckets[suffixes[i].bid].size++;
        ++bucket_size_total;
    }
    /* last iteration */
    {
        size_t i=total_amount_to_recv-1;
        assert(suffixes[i].sid < last_id);
        assert(suffixes[i].bid < n_buckets);
        assert(suffixes[i].next == NULL);
        suffixes[i].next = NULL;
        if (buckets[suffixes[i].bid].suffixes == NULL) {
            buckets[suffixes[i].bid].suffixes = &suffixes[i];
            bucket_offset[suffixes[i].bid] = i;
        }
        buckets[suffixes[i].bid].size++;
        ++bucket_size_total;
    }
    mpix_allreduce(bucket_offset, MPI_SUM, comm);
#else
    /* The suffixes contains sorted suffixes based on the buckets they belong
     * to. That means we can simply update the 'next' links! */
    size_t last_id = sequences->get_global_count();
    bucket_size_total = 0;
    for (size_t i=0; i<total_amount_to_recv; ++i) {
        assert(suffixes[i].sid < last_id);
        assert(suffixes[i].bid < n_buckets);
        assert(suffixes[i].next == NULL);
        suffixes[i].next = buckets[suffixes[i].bid].suffixes;
        buckets[suffixes[i].bid].suffixes = &suffixes[i];
        buckets[suffixes[i].bid].size++;
        ++bucket_size_total;
    }
#endif

    suffixes_size = n_suffixes;
    buckets_size = n_buckets;

#if 0
    bool found_start = false;
    bool found_stop = false;
    for (size_t i=0; i<buckets_size; ++i) {
        if (found_start && found_stop) {
            assert(bucket_owner[i] != comm_rank);
        }
        else if (found_start) {
            if (bucket_owner[i] != comm_rank) {
                found_stop = true;
                last_bucket = i-1;
            }
        }
        else if (found_stop) {
            assert(0);
        }
        else {
            if (bucket_owner[i] == comm_rank) {
                found_start = true;
                first_bucket = i;
            }
        }
    }
    if (found_start && !found_stop) {
        last_bucket = buckets_size-1;
    }
#endif

#if 1
    mpix_print_sync("first_bucket", first_bucket, comm);
    mpix_print_sync("last_bucket", last_bucket, comm);
    mpix_print_sync("bucket_size_total", bucket_size_total, comm);
#endif
}


SuffixBuckets::~SuffixBuckets()
{
#if USE_ARMCI
    ARMCI_Free(suffixes);
#else
    delete [] suffixes;
#endif
    delete [] buckets;
}


bool SuffixBuckets::owns(size_t bid) const
{
    assert(bid < buckets_size);

    size_t owner = bucket_owner[bid];

    return owner == comm_rank;
}


Suffix* SuffixBuckets::get(size_t bid)
{
#if USE_ARMCI
    assert(bid < buckets_size);

    size_t size = bucket_size[bid];
    size_t owner = bucket_owner[bid];
    size_t offset = bucket_offset[bid];
    void *address = bucket_address[owner];

    if (owner == comm_rank) {
        /* already owned, just return it */
        return buckets[bid].suffixes;
    }
    else {
        LockGuard<PthreadMutex> guard(mutex);
        //Suffix *remote_suffixes = new Suffix[size];
        Suffix *remote_suffixes = (Suffix*)ARMCI_Malloc_local(sizeof(Suffix)*size);
        ARMCI_Get((void*)&bucket_address[owner][offset],
                remote_suffixes, size*sizeof(Suffix), owner);
        for (size_t i=0; i<size-1; ++i) {
            assert(remote_suffixes[i].bid == bid);
            remote_suffixes[i].next = &remote_suffixes[i+1];
        }
        assert(remote_suffixes[size-1].next == NULL);
        return remote_suffixes;
    }
#else
    assert(0);
#endif
}


}; /* namespace pgraph */

