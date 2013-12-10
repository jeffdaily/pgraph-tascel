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
#include <vector>

#include <tascel.h>

#include "constants.h"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"

#if HAVE_ARMCI
#define USE_ARMCI 1
#else
#define USE_ARMCI 0
#endif

#if USE_ARMCI
/* ANL's armci does not extern "C" inside the header */
extern "C" {
#include <armci.h>
}
#endif

using tascel::LockGuard;
using tascel::PthreadMutex;
using std::vector;

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
        const char tmp = kmer[i] - 'A';
        if (tmp < 0) {
            printf("kmer[%d]=(%c) - 'A' < 0\n", i, kmer[i]);
        }
        if (tmp >= SIGMA) {
            printf("kmer[%d]=(%c) >= SIGMA=(%u)\n", i, kmer[i], SIGMA);
        }
        assert(tmp >= 0);
        assert(tmp < SIGMA);
        value = value * SIGMA + tmp;
    }

    return value;
}


struct SuffixBucketIndexFunctor {
    SuffixBucketIndexFunctor(vector<size_t> *bucket_owner)
        :   bucket_owner(bucket_owner)
    { }

    bool operator()(const Suffix &i, const Suffix &j) {
        return (*bucket_owner)[i.bid] < (*bucket_owner)[j.bid];
    }

    vector<size_t> *bucket_owner;
};


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
                             const SplitEnum &split_type)
    :   comm_rank(-1)
    ,   comm_size(-1)
    ,   sequences(sequences)
    ,   param(param)
    ,   suffixes(NULL)
    ,   suffixes_size(0)
    ,   buckets(NULL)
    ,   buckets_size(0)
    ,   my_buckets()
    ,   bucket_size()
    ,   bucket_owner()
    ,   bucket_offset()
    ,   bucket_address(NULL)
    ,   bucket_size_total(0)
    ,   mutex()
    ,   split_type(split_type)
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
    vector<pair<size_t,size_t> > bucket_sorted_map(n_buckets);

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
#if DEBUG || 1
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
            if (bucket_index >= n_buckets) {
                printf("[%d] bucket_index >= n_buckets (%zu >= %zu)\n",
                        comm_rank, bucket_index, n_buckets);
                printf("[%d] j=%zu stop_index=%zu k=%d data='%s' len=%zu\n",
                        comm_rank, j, stop_index, param.window_size,
                        std::string(sequence_data,sequence_length).c_str(),
                        sequence_length);
            }
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

    if (split_type == SPLIT_DUMB) {
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
    else if (split_type == SPLIT_SUFFIXES) {
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
            }
            bucket_owner[i] = rank;
        }
    }
    else if (split_type == SPLIT_SORTED) {
        for (size_t i=0; i<n_buckets; ++i) {
            bucket_sorted_map[i].first = bucket_size[i];
            bucket_sorted_map[i].second = i;
        }
        std::sort(bucket_sorted_map.begin(), bucket_sorted_map.end());
        int rank = 0;
        for (size_t i=0; i<n_buckets; ++i) {
            bucket_owner[bucket_sorted_map[i].second] = rank;
            rank = (rank + 1) % comm_size;
        }
    }
    else {
        assert(0);
    }
#if DEBUG
    mpix_print_zero("bucket_owner", vec_to_string(bucket_owner), comm);
#endif
    for (size_t i=0; i<n_buckets; ++i) {
        if (bucket_owner[i] == comm_rank) {
            my_buckets.push_back(i);
        }
    }
#if DEBUG
    mpix_print_sync("my_buckets", vec_to_string(my_buckets), comm);
#endif

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
    mpix_print_sync("total_amount_to_recv", total_amount_to_recv, comm);

    /* We are preparing for the all to all, so we sort the suffixes. Yes, this
     * invalidates any of their 'next' pointers but these will be rebuilt after
     * the all to all. */
    std::sort(initial_suffixes.begin(),
              initial_suffixes.end(),
              //SuffixBucketIndexCompare
              SuffixBucketIndexFunctor(&bucket_owner)
    );
#if USE_ARMCI
    bucket_address = new Suffix*[comm_size];
    (void)ARMCI_Malloc((void**)bucket_address,
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
        MPI_UNSIGNED_LONG
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

    if (total_amount_to_recv > 0) {
        std::sort(suffixes,
                suffixes+total_amount_to_recv,
                SuffixBucketIndexCompare);
        /* The suffixes contains sorted suffixes based on the buckets they belong
         * to. That means we can simply update the 'next' links! */
        size_t last_id = sequences->get_global_count();
        bucket_size_total = 0;
        for (size_t i=0,limit=total_amount_to_recv-1; i<limit; ++i) {
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
    }
    mpix_allreduce(bucket_offset, MPI_SUM, comm);

    suffixes_size = n_suffixes;
    buckets_size = n_buckets;

#if 1
    mpix_print_sync("bucket_size_total", bucket_size_total, comm);
#endif
}


SuffixBuckets::~SuffixBuckets()
{
#if USE_ARMCI
    ARMCI_Free(suffixes);
    delete [] bucket_address;
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
    assert(owns(bid));
    return buckets[bid].suffixes;
#endif
}


#if USE_ARMCI

struct SuffixBucket2IndexFunctor {
    SuffixBucket2IndexFunctor(size_t comm_size)
        :   comm_size(comm_size)
    { }

    bool operator()(const Suffix &i, const Suffix &j) {
        size_t owner_i = i.bid % comm_size;
        size_t owner_j = j.bid % comm_size;
        return owner_i < owner_j;
    }

    size_t comm_size;
};

#define DEBUG 1

SuffixBuckets2::SuffixBuckets2(SequenceDatabase *sequences,
                             const Parameters &param,
                             MPI_Comm comm)
    :   comm_rank(-1)
    ,   comm_size(-1)
    ,   n_buckets(0)
    ,   sequences(sequences)
    ,   param(param)
    ,   suffixes_remote(NULL)
    ,   suffixes(NULL)
    ,   suffixes_size(0)
    ,   buckets_remote(NULL)
    ,   buckets(NULL)
    ,   buckets_size()
    ,   mutex()
    ,   count_remote_buckets(0)
    ,   count_remote_suffixes(0)
{
    size_t n_suffixes = 0;
    size_t suffix_index = 0;
    size_t n_seq = 0;
    size_t remainder = 0;
    size_t start = 0;
    size_t stop = 0;
    int ierr = 0;

    assert(USE_ARMCI);

    ierr = MPI_Comm_rank(comm, &comm_rank);
    assert(MPI_SUCCESS == ierr);
    ierr = MPI_Comm_size(comm, &comm_size);
    assert(MPI_SUCCESS == ierr);

    /* allocate buckets */
    n_buckets = powz(SIGMA, param.window_size);
    mpix_print_zero("n_buckets", n_buckets, comm);

    /* how many buckets does this process own? */
    buckets_size = n_buckets / comm_size;
    if (comm_rank < (n_buckets % comm_size)) {
        buckets_size += 1;
    }

    /* allocate suffixes */
    n_suffixes = sequences->get_global_size()
        - sequences->get_global_count() * param.window_size;
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
        initial_suffixes_size +=
            sequence.get_sequence_length() - param.window_size;
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
            suffix_index++;
        }
    }
    initial_suffixes.resize(suffix_index);

#if DEBUG
    mpix_print_sync("suffix_index", suffix_index, comm);
#endif

    // Note: bucket owner is calculated as bucket ID % comm_size

    vector<int> amount_to_send(comm_size, 0);
    vector<int> amount_to_recv(comm_size, 0);
    for (size_t i=0; i<initial_suffixes.size(); ++i) {
        size_t owner = initial_suffixes[i].bid % comm_size;
        amount_to_send[owner] += 1;
    }
#if DEBUG
    size_t desitinations = 0;
    for (size_t i=0; i<comm_size; ++i) {
        if (amount_to_send[i] > 0) {
            desitinations += 1;
        }
    }
    mpix_print_sync("desitinations", desitinations, comm);
#endif
#if DEBUG
    mpix_print_sync("amount_to_send", vec_to_string(amount_to_send), comm);
#endif
    mpix_alltoall(amount_to_send, amount_to_recv, comm);
#if DEBUG
    mpix_print_sync("amount_to_recv", vec_to_string(amount_to_recv), comm);
#endif

    int total_amount_to_send = 0;
    int total_amount_to_recv = 0;
    for (int i=0; i<comm_size; ++i) {
        total_amount_to_send += amount_to_send[i];
        total_amount_to_recv += amount_to_recv[i];
    }
    mpix_print_sync("total_amount_to_send", total_amount_to_send, comm);
    mpix_print_sync("total_amount_to_recv", total_amount_to_recv, comm);
    suffixes_size = total_amount_to_recv;

    /* We are preparing for the all to all, so we sort the suffixes.
     * Sort is based on owner ID. */
    std::sort(initial_suffixes.begin(),
              initial_suffixes.end(),
              SuffixBucket2IndexFunctor(comm_size)
    );

    suffixes_remote = new Suffix*[comm_size];
    (void)ARMCI_Malloc((void**)suffixes_remote,
            sizeof(Suffix)*suffixes_size);
    suffixes = suffixes_remote[comm_rank];
    mpix_print_sync("suffixes_size", suffixes_size, comm);

    buckets_remote = new Bucket2*[comm_size];
    (void)ARMCI_Malloc((void**)buckets_remote,
            sizeof(Bucket2)*buckets_size);
    buckets = buckets_remote[comm_rank];
    mpix_print_sync("buckets_size", buckets_size, comm);
    for (size_t i=0; i<buckets_size; ++i) {
        buckets[i].offset = 0;
        buckets[i].size = 0;
    }

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

    /* need to alltoallv the buckets to the owning processes */
    MPI_Datatype SuffixType;
    MPI_Datatype type[4] = {
        mpix_get_mpi_datatype(initial_suffixes[0].sid),
        mpix_get_mpi_datatype(initial_suffixes[0].pid),
        mpix_get_mpi_datatype(initial_suffixes[0].bid),
        MPI_UNSIGNED_LONG /* void* */
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
#if DEBUG
    mpix_print_sync("initial_suffixes", vec_to_string(initial_suffixes), comm);
#endif
    ierr = MPI_Alltoallv(
            &initial_suffixes[0], &amount_to_send[0],
            &send_displacements[0], SuffixType,
            &suffixes[0], &amount_to_recv[0],
            &recv_displacements[0], SuffixType, comm);
    assert(MPI_SUCCESS == ierr);
#if DEBUG
    mpix_print_sync("suffixes", arr_to_string(suffixes, total_amount_to_recv), comm);
#endif

    size_t bucket_size_total = 0;
    if (total_amount_to_recv > 0) {
        std::sort(suffixes,
                suffixes+total_amount_to_recv,
                SuffixBucketIndexCompare);
        /* The suffixes contains sorted suffixes based on the buckets they belong
         * to. That means we can simply update the 'next' links! */
        size_t last_id = sequences->get_global_count();
        for (size_t i=0,limit=total_amount_to_recv-1; i<limit; ++i) {
            assert(suffixes[i].sid < last_id);
            assert(suffixes[i].bid < n_buckets);
            assert(suffixes[i].bid <= suffixes[i+1].bid);
            assert(suffixes[i].bid % comm_size == comm_rank);
            assert(suffixes[i].next == NULL);
            if (suffixes[i].bid == suffixes[i+1].bid) {
                suffixes[i].next = &suffixes[i+1];
            }
            else {
                suffixes[i].next = NULL;
            }
            size_t bucket_index = suffixes[i].bid / comm_size;
            if (buckets[bucket_index].size == 0) {
                buckets[bucket_index].offset = i;
            }
            buckets[bucket_index].size++;
            ++bucket_size_total;
        }
        /* last iteration */
        {
            size_t i=total_amount_to_recv-1;
            assert(suffixes[i].sid < last_id);
            assert(suffixes[i].bid < n_buckets);
            assert(suffixes[i].bid % comm_size == comm_rank);
            assert(suffixes[i].next == NULL);
            suffixes[i].next = NULL;
            size_t bucket_index = suffixes[i].bid / comm_size;
            if (buckets[bucket_index].size == 0) {
                buckets[bucket_index].offset = i;
            }
            buckets[bucket_index].size++;
            ++bucket_size_total;
        }
    }

#if 1
    mpix_print_sync("bucket_size_total", bucket_size_total, comm);
#endif
}


SuffixBuckets2::~SuffixBuckets2()
{
    ARMCI_Free(suffixes);
    ARMCI_Free(buckets);
    delete [] suffixes_remote;
    delete [] buckets_remote;
}


bool SuffixBuckets2::owns(size_t bid) const
{
    assert(bid < n_buckets);

    return (bid % comm_size) == comm_rank;
}


Bucket* SuffixBuckets2::get(size_t bid)
{
    size_t owner = bid % comm_size;
    size_t bucket_index = bid / comm_size;
    Bucket *bucket = new Bucket;

    if (owner == comm_rank) {
        /* already owned, just return it */
        bucket->suffixes = &suffixes[buckets[bucket_index].offset];
        bucket->size = buckets[bucket_index].size;
    }
    else {
        LockGuard<PthreadMutex> guard(mutex);
        /* get remote bucket info */
        Bucket2 *remote_bucket = (Bucket2*)ARMCI_Malloc_local(sizeof(Bucket2));
        ARMCI_Get(&buckets_remote[owner][bucket_index],
                remote_bucket, sizeof(Bucket2), owner);
        bucket->size = remote_bucket->size;
        count_remote_buckets += 1;
        count_remote_suffixes += bucket->size;
        if (bucket->size > 0) {
            /* get remote suffixes */
            Suffix *remote_suffixes =
                (Suffix*)ARMCI_Malloc_local(
                        sizeof(Suffix)*remote_bucket->size);
            ARMCI_Get(&suffixes_remote[owner][remote_bucket->offset],
                    remote_suffixes, sizeof(Suffix)*remote_bucket->size, owner);
            for (size_t i=0; i<remote_bucket->size-1; ++i) {
                assert(remote_suffixes[i].bid == bid);
                remote_suffixes[i].next = &remote_suffixes[i+1];
            }
            assert(remote_suffixes[remote_bucket->size-1].next == NULL);
            bucket->suffixes = remote_suffixes;
        }
        else {
            bucket->suffixes = NULL;
        }
        ARMCI_Free_local(remote_bucket);
    }

    return bucket;
}


void SuffixBuckets2::rem(size_t bid, Bucket *bucket)
{
    if (!owns(bid) && bucket->size > 0 && bucket->suffixes != NULL) {
        LockGuard<PthreadMutex> guard(mutex);
        ARMCI_Free_local(bucket->suffixes);
    }
    delete bucket;
}

#endif /* USE_ARMCI */


}; /* namespace pgraph */

