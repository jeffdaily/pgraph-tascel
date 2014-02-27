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
#include <vector>

#include <tascel.h>

#include "constants.h"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixBucketsTascel.hpp"

using namespace tascel;
using ::std::vector;

extern Dispatcher<NullMutex> serverDispatcher;

namespace pgraph {

SuffixBucketsTascel::SuffixBucketsTascel(SequenceDatabase *sequences,
                             const Parameters &param,
                             MPI_Comm comm)
    :   SuffixBuckets(sequences, param, comm)
    ,   aid_suffixes()
    ,   aid_meta()
    ,   suffixes(NULL)
    ,   suffixes_size(0)
    ,   buckets(NULL)
    ,   buckets_size(0)
    ,   n_nonempty(0)
    ,   count_remote_buckets(0)
    ,   count_remote_suffixes(0)
    ,   owned_buckets()
{
    size_t n_suffixes = 0;
    size_t suffix_index = 0;
    size_t n_seq = 0;
    size_t remainder = 0;
    size_t start = 0;
    size_t stop = 0;
    int ierr = 0;

    /* how many buckets does this process own? */
    buckets_size = n_buckets / comm_size;
    if (comm_rank < (n_buckets % comm_size)) {
        buckets_size += 1;
    }

    /* how many suffixes are in the given sequence database? */
    n_suffixes = sequences->char_size()
        - sequences->size() * param.window_size;

#if DEBUG || 1
    mpix_print_zero("n_buckets", n_buckets, comm);
    mpix_print_sync("buckets_size", buckets_size, comm);
    mpix_print_zero("n_suffixes", n_suffixes, comm);
#endif

    /* each MPI rank gets a contiguous range of sequences to bucket */
    n_seq = sequences->size() / comm_size;
    remainder = sequences->size() % comm_size;
    start = n_seq * comm_rank;
    stop = n_seq * (comm_rank+1);
    if (comm_rank < remainder) {
        start += comm_rank;
        stop += comm_rank+1;
        n_seq += 1;
    }
    else {
        start += remainder;
        stop += remainder;
    }
    if (stop > sequences->size()) {
        stop = sequences->size();
    }
#if DEBUG || 1
    mpix_print_sync("n_seq", n_seq, comm);
    mpix_print_sync("remainder", remainder, comm);
    mpix_print_sync("start", start, comm);
    mpix_print_sync("stop", stop, comm);
#endif

    size_t initial_suffixes_size = 0;
    for (size_t i = start; i < stop; ++i) {
        initial_suffixes_size += (*sequences)[i].size() - param.window_size;
    }
    vector<Suffix> initial_suffixes(initial_suffixes_size);

    /* slide k-mers for every sequence and bucket them */
    for (size_t i = start; i < stop; ++i) {
        size_t stop_index = 0;
        const char *sequence_data = NULL;
        size_t sequence_length = 0;

        (*sequences)[i].get_sequence(sequence_data, sequence_length);
        stop_index = sequence_length - param.window_size - 1;

        if (sequence_length <= ((unsigned)param.window_size)) continue;

        for (size_t j = 0; j <= stop_index; ++j) {
            if (filter_out(sequence_data + j)) {
                printf("[%zu] filtered out '%s'\n",
                        comm_rank, string(sequence_data+j,param.window_size).c_str());
            } else {
                size_t bid = bucket_index(sequence_data + j);
                assert(bid < n_buckets);
                initial_suffixes[suffix_index].sid = i;
                initial_suffixes[suffix_index].pid = j;
                initial_suffixes[suffix_index].bid = bid;
                initial_suffixes[suffix_index].next = NULL;
                suffix_index++;
            }
        }
    }
    initial_suffixes.resize(suffix_index);

#if DEBUG || 1
    mpix_print_sync("suffix_index", suffix_index, comm);
#endif

    // Note: bucket owner is calculated as bucket ID % comm_size

    vector<int> amount_to_send(comm_size, 0);
    vector<int> amount_to_recv(comm_size, 0);
    for (size_t i=0; i<initial_suffixes.size(); ++i) {
        size_t owner = initial_suffixes[i].bid % comm_size;
        amount_to_send[owner] += 1;
    }
#if DEBUG || 1
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
    for (size_t i=0; i<comm_size; ++i) {
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
              SuffixOwnerCompareFunctor(comm_size)
    );

    aid_suffixes = theRma().allocColl(sizeof(Suffix)*suffixes_size);
    suffixes = reinterpret_cast<Suffix*>(
            theRma().lookupPointer(RmaPtr(aid_suffixes)));
    mpix_print_sync("suffixes_size", suffixes_size, comm);

    aid_meta = theRma().allocColl(sizeof(BucketMeta)*buckets_size);
    buckets = reinterpret_cast<BucketMeta*>(
            theRma().lookupPointer(RmaPtr(aid_meta)));
    mpix_print_sync("buckets_size", buckets_size, comm);
    for (size_t i=0; i<buckets_size; ++i) {
        size_t bid = i*comm_size + comm_rank;
        buckets[i].offset = 0;
        buckets[i].size = 0;
        buckets[i].bid = bid;
        owned_buckets.push_back(bid);
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
                SuffixBucketCompare);
        /* The suffixes contains sorted suffixes based on the buckets they belong
         * to. That means we can simply update the 'next' links! */
        size_t last_id = sequences->size();
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
                n_nonempty++;
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
                n_nonempty++;
            }
            buckets[bucket_index].size++;
            ++bucket_size_total;
        }
    }

    mpix_allreduce(n_nonempty, MPI_SUM, comm);
#if 1
    mpix_print_sync("bucket_size_total", bucket_size_total, comm);
#endif
}


SuffixBucketsTascel::~SuffixBucketsTascel()
{
    theRma().deallocColl(aid_suffixes);
    theRma().deallocColl(aid_meta);
}


bool SuffixBucketsTascel::owns(size_t bid) const
{
    assert(bid < n_buckets);

    return (bid % comm_size) == comm_rank;
}


Bucket* SuffixBucketsTascel::get(size_t bid)
{
    size_t owner = bid % comm_size;
    size_t bucket_index = bid / comm_size;
    Bucket *bucket = new Bucket;

    if (owner == comm_rank) {
        /* already owned, just return it */
        bucket->size = buckets[bucket_index].size;
        bucket->bid = bid;
        if (bucket->size > 0) {
            bucket->suffixes = &suffixes[buckets[bucket_index].offset];
        }
        else {
            bucket->suffixes = NULL;
        }
    }
    else {
        /* get remote bucket info */
        RmaRequest *localReq = RmaRequest::construct();
        RmaRequest *remoteReq = RmaRequest::construct();
        BucketMeta *remote_bucket = new BucketMeta;
        theRma().get(remote_bucket,
                RmaPtr(aid_meta, bucket_index*sizeof(BucketMeta)),
                sizeof(BucketMeta),
                owner, localReq, remoteReq);
        while(!remoteReq->test()) {
            AmListenObjCodelet<NullMutex>* lcodelet;
            if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                lcodelet->execute();
            }
        }
        delete remoteReq;
        while(!localReq->test()) {
            AmListenObjCodelet<NullMutex>* lcodelet;
            if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                lcodelet->execute();
            }
        }
        delete localReq;
        assert(remote_bucket->bid == bid);

        bucket->size = remote_bucket->size;
        count_remote_buckets += 1;
        count_remote_suffixes += bucket->size;
        if (bucket->size > 0) {
            /* get remote suffixes */
            localReq = RmaRequest::construct();
            remoteReq = RmaRequest::construct();
            Suffix *remote_suffixes = new Suffix[remote_bucket->size];
            theRma().get(remote_suffixes,
                    RmaPtr(aid_suffixes, remote_bucket->offset),
                    sizeof(Suffix)*remote_bucket->size,
                    owner, localReq, remoteReq);
            while(!remoteReq->test()) {
                AmListenObjCodelet<NullMutex>* lcodelet;
                if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                    lcodelet->execute();
                }
            }
            delete remoteReq;
            while(!localReq->test()) {
                AmListenObjCodelet<NullMutex>* lcodelet;
                if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                    lcodelet->execute();
                }
            }
            delete localReq;

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
        bucket->bid = bid;
        delete remote_bucket;
    }

    return bucket;
}


void SuffixBucketsTascel::rem(Bucket *bucket)
{
    if (!owns(bucket->bid) && bucket->size > 0 && bucket->suffixes != NULL) {
        delete [] bucket->suffixes;
    }
    delete bucket;
}

}; /* namespace pgraph */

