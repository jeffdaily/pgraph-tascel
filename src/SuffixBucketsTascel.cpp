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

#include "mpix.hpp"
#include "mpix-types.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixBucketsTascel.hpp"

using namespace tascel;
using ::std::vector;

extern Dispatcher<NullMutex> serverDispatcher;

namespace pgraph {

/** Functor for sorting Suffix instances based on Bucket owner. */
struct SuffixOwnerCompareFunctor {
    size_t comm_size;

    SuffixOwnerCompareFunctor(size_t comm_size)
        :   comm_size(comm_size)
    { }

    bool operator()(const Suffix &i, const Suffix &j) {
        size_t owner_i = i.bid % comm_size;
        size_t owner_j = j.bid % comm_size;
        return owner_i < owner_j;
    }
};

static bool SuffixBucketCompare(const Suffix &i, const Suffix &j) {
    return i.bid < j.bid;
}


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
{
    size_t n_buckets = powz(SIGMA,param.window_size);
    size_t n_suffixes = 0;
    size_t n_suffixes_skipped = 0;
    size_t suffix_index = 0;
    size_t n_seq = 0;
    size_t remainder = 0;
    size_t start = 0;
    size_t stop = 0;
    int ierr = 0;

    /* how many suffixes are in the given sequence database? */
    n_suffixes = sequences->char_size()
        - sequences->size() * param.window_size;

#if DEBUG || 1
    mpix::print_zero("n_buckets", n_buckets, comm);
#endif
#if DEBUG || 1
    mpix::print_zero("n_suffixes", n_suffixes, comm);
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
#if DEBUG
    mpix::print_sync("n_seq", n_seq, comm);
    mpix::print_sync("remainder", remainder, comm);
    mpix::print_sync("start", start, comm);
    mpix::print_sync("stop", stop, comm);
#endif

    size_t initial_suffixes_size = 0;
    for (size_t i = start; i < stop; ++i) {
        size_t sequence_length = (*sequences)[i].size();
        if (sequence_length >= param.window_size) {
            initial_suffixes_size += sequence_length - param.window_size + 1;
        }
    }
    vector<Suffix> initial_suffixes(initial_suffixes_size);

    /* slide k-mers for every sequence and bucket them */
    for (size_t i = start; i < stop; ++i) {
        size_t stop_index = 0;
        const char *sequence_data = NULL;
        size_t sequence_length = 0;

        (*sequences)[i].get_sequence(sequence_data, sequence_length);
        if (sequence_length < ((unsigned)param.window_size)) {
            n_suffixes_skipped += 1;
            continue;
        }
        stop_index = sequence_length - param.window_size;
        for (size_t j = 0; j <= stop_index; ++j) {
            if (filter_out(sequence_data + j)) {
                n_suffixes_skipped += 1;
#if 0
                printf("[%zu] filtered out '%s'\n",
                        comm_rank, string(sequence_data+j,param.window_size).c_str());
#endif
            } else {
                size_t bid = bucket_index(sequence_data + j);
                assert(bid < n_buckets);
                initial_suffixes[suffix_index].sid = i;
                initial_suffixes[suffix_index].pid = j;
                initial_suffixes[suffix_index].bid = bid;
                initial_suffixes[suffix_index].k = param.window_size;
                initial_suffixes[suffix_index].next = NULL;
                suffix_index++;
            }
        }
    }
    initial_suffixes.resize(suffix_index);
#if DEBUG || 1
    mpix::print_sync("n_suffixes_skipped", n_suffixes_skipped, comm);
#endif

#if DEBUG
    mpix::print_sync("suffix_index", suffix_index, comm);
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
    mpix::print_sync("desitinations", desitinations, comm);
#endif
#if DEBUG
    mpix::print_sync("amount_to_send", amount_to_send, comm);
#endif

    mpix::alltoall(amount_to_send, amount_to_recv, comm);

#if DEBUG
    mpix::print_sync("amount_to_recv", amount_to_recv, comm);
#endif

    int total_amount_to_send = 0;
    int total_amount_to_recv = 0;
    for (size_t i=0; i<comm_size; ++i) {
        total_amount_to_send += amount_to_send[i];
        total_amount_to_recv += amount_to_recv[i];
    }
#if DEBUG
    mpix::print_sync("total_amount_to_send", total_amount_to_send, comm);
    mpix::print_sync("total_amount_to_recv", total_amount_to_recv, comm);
#endif
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
#if DEBUG
    mpix::print_sync("suffixes_size", suffixes_size, comm);
#endif

    vector<int> send_displacements(comm_size, 0);
    vector<int> recv_displacements(comm_size, 0);

    for (size_t i=1; i<comm_size; ++i) {
        send_displacements[i] = send_displacements[i-1] + amount_to_send[i-1];
        recv_displacements[i] = recv_displacements[i-1] + amount_to_recv[i-1];
    }
#if DEBUG
    mpix::print_sync("send_displacements", send_displacements, comm);
    mpix::print_sync("recv_displacements", recv_displacements, comm);
#endif

    /* need to alltoallv the buckets to the owning processes */
    MPI_Datatype SuffixType = mpix::get_mpi_datatype<Suffix>(initial_suffixes[0]);
#if DEBUG
    mpix::print_sync("initial_suffixes", initial_suffixes, comm);
#endif

    ierr = MPI_Alltoallv(
            &initial_suffixes[0], &amount_to_send[0],
            &send_displacements[0], SuffixType,
            &suffixes[0], &amount_to_recv[0],
            &recv_displacements[0], SuffixType, comm);
    assert(MPI_SUCCESS == ierr);

#if DEBUG
    mpix::print_sync("suffixes", suffixes, total_amount_to_recv, comm);
#endif

    vector<BucketMeta> initial_buckets;

    if (total_amount_to_recv > 0) {
        /* sort the received Suffix instances based on bucket ID */
        std::sort(suffixes, suffixes+total_amount_to_recv,
                SuffixBucketCompare);
        /* The suffixes contains sorted suffixes based on the buckets
         * they belong to. That means we can simply update the 'next'
         * links! */
        size_t last_id = sequences->size();
        size_t longest = sequences->longest();
        for (size_t i=0; i<total_amount_to_recv; ++i) {
            assert(suffixes[i].sid < last_id);
            assert(suffixes[i].pid < longest);
            assert(suffixes[i].bid < n_buckets);
            assert(suffixes[i].bid % comm_size == comm_rank);
            assert(suffixes[i].k == param.window_size);
            assert(suffixes[i].next == NULL);
            if (i == 0 || suffixes[i-1].bid != suffixes[i].bid) {
                initial_buckets.push_back(BucketMeta(
                            i, 1, suffixes[i].bid, suffixes[i].k));
                n_nonempty++;
            }
            else if (suffixes[i-1].bid == suffixes[i].bid) {
                suffixes[i-1].next = &suffixes[i];
                initial_buckets.back().size += 1;
            }
            else {
                assert(0); /* this shouldn't be reached */
            }
        }
    }

    /* now that initial bucket are created, collect sizes stats */
    for (size_t i=0; i<initial_buckets.size(); ++i) {
        buckets_stats.push_back(initial_buckets[i].size);
    }
    {
        Stats *all_stats = new Stats[comm_size];
        mpix::check(MPI_Gather(&buckets_stats, sizeof(Stats), MPI_CHAR,
                    all_stats, sizeof(Stats), MPI_CHAR,
                    0, comm));
        if (0 == comm_rank) {
            Stats combined_stats;
            cout << Stats::header() << endl;
            for (size_t i=0; i<comm_size; ++i) {
                combined_stats.push_back(all_stats[i]);
                cout << all_stats[i] << endl;
            }
            buckets_stats = combined_stats;
        }
        mpix::bcast(buckets_stats, 0, comm);
        delete [] all_stats;
    }

#if DEBUG || 1
    mpix::print_zero(" stats header", Stats::header(), comm);
    mpix::print_zero("buckets_stats", buckets_stats, comm);
#endif

#if DEBUG
    {
        Stats *all_stats = new Stats[comm_size];
        mpix::check(MPI_Gather(&buckets_stats, sizeof(Stats), MPI_CHAR,
                    all_stats, sizeof(Stats), MPI_CHAR,
                    0, comm));
        if (0 == comm_rank) {
            cout << Stats::header() << endl;
            for (size_t i=0; i<comm_size; ++i) {
                cout << all_stats[i] << endl;
            }
        }
        delete [] all_stats;
    }
#endif

    /* now that initial buckets are created, refine the big ones */
    bool done = false;
    while (!done) {
        done = true;
        for (size_t i=0,limit=initial_buckets.size(); i<limit; ++i) {
            BucketMeta &bucket = initial_buckets[i];
            if (bucket.size > 2*buckets_stats.stddev()
                    && bucket.k < param.window_size) {
                refine_bucket(bucket, initial_buckets);
                done = false;
            }
        }
    }
    buckets_size = initial_buckets.size();

#if DEBUG || 1
    mpix::print_sync("buckets_size", buckets_size, comm);
#endif
    aid_meta = theRma().allocColl(sizeof(BucketMeta)*buckets_size);
    buckets = reinterpret_cast<BucketMeta*>(
            theRma().lookupPointer(RmaPtr(aid_meta)));
    ::std::copy(initial_buckets.begin(), initial_buckets.end(), buckets);

    buckets_size_global = buckets_size;
    mpix::allreduce(buckets_size_global, MPI_SUM, comm);
#if DEBUG || 1
    mpix::print_sync("buckets_size_global", buckets_size_global, comm);
#endif
    mpix::allreduce(n_nonempty, MPI_SUM, comm);
#if DEBUG || 1
    mpix::print_sync("n_nonempty", n_nonempty, comm);
#endif
}


SuffixBucketsTascel::~SuffixBucketsTascel()
{
    theRma().deallocColl(aid_suffixes);
    theRma().deallocColl(aid_meta);
}


Bucket* SuffixBucketsTascel::get(int owner, size_t bucket_index)
{
    Bucket *bucket = new Bucket;

    if (owner == comm_rank) {
        /* already owned, just return it */
        bucket->size = buckets[bucket_index].size;
        bucket->bid = buckets[bucket_index].bid;
        bucket->k = bucket[bucket_index].k;
        bucket->owner = owner;
        if (bucket->size > 0) {
            bucket->suffixes = &suffixes[buckets[bucket_index].offset];
        }
        else {
            bucket->suffixes = NULL;
        }
    }
    else {
        /* get remote bucket info */
        Dispatcher<NullMutex> dispatcher;
        RmaRequest *localReq = RmaRequest::construct();
        RmaRequest *remoteReq = RmaRequest::construct();
        BucketMeta *remote_bucket = new BucketMeta;
        theRma().get(remote_bucket,
                RmaPtr(aid_meta, bucket_index*sizeof(BucketMeta)),
                sizeof(BucketMeta),
                owner, localReq, remoteReq);
        dispatcher.registerCodelet(localReq);
        dispatcher.registerCodelet(remoteReq);
        while (!dispatcher.empty()) {
            Codelet* codelet;
            if ((codelet = dispatcher.progress()) != NULL) {
                codelet->execute();
            }
#if !defined(THREADED)
            AmListenObjCodelet<NullMutex>* lcodelet;
            if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                lcodelet->execute();
            }
#endif
        }
        delete remoteReq;
        delete localReq;

        bucket->size = remote_bucket->size;
        bucket->bid = remote_bucket->bid;
        bucket->k = remote_bucket->k;
        bucket->owner = owner;
        count_remote_buckets += 1;
        count_remote_suffixes += bucket->size;
        if (bucket->size > 0) {
            /* get remote suffixes */
            localReq = RmaRequest::construct();
            remoteReq = RmaRequest::construct();
            Suffix *remote_suffixes = new Suffix[remote_bucket->size];
            theRma().get(remote_suffixes,
                    RmaPtr(aid_suffixes, sizeof(Suffix)*remote_bucket->offset),
                    sizeof(Suffix)*remote_bucket->size,
                    owner, localReq, remoteReq);
            dispatcher.registerCodelet(localReq);
            dispatcher.registerCodelet(remoteReq);
            while (!dispatcher.empty()) {
                Codelet* codelet;
                if ((codelet = dispatcher.progress()) != NULL) {
                    codelet->execute();
                }
#if !defined(THREADED)
                AmListenObjCodelet<NullMutex>* lcodelet;
                if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                    lcodelet->execute();
                }
#endif
            }
            delete remoteReq;
            delete localReq;

            for (size_t i=0; i<remote_bucket->size-1; ++i) {
                assert(remote_suffixes[i].bid == remote_suffixes[i+1].bid);
                remote_suffixes[i].next = &remote_suffixes[i+1];
            }
            assert(remote_suffixes[remote_bucket->size-1].next == NULL);
            bucket->suffixes = remote_suffixes;
        }
        else {
            bucket->suffixes = NULL;
        }
        delete remote_bucket;
    }

    return bucket;
}


void SuffixBucketsTascel::rem(Bucket *bucket)
{
    if (bucket->owner != comm_rank && bucket->size > 0 && bucket->suffixes != NULL) {
        delete [] bucket->suffixes;
    }
    delete bucket;
}


void SuffixBucketsTascel::refine_bucket(
        BucketMeta &bucket,
        vector<BucketMeta> &buckets)
{
    Suffix *local_suffixes = &suffixes[bucket.offset];
    size_t local_suffixes_size = bucket.size;

    for (size_t i=0; i<local_suffixes_size; ++i) {
        const char *sequence_data = NULL;
        size_t sequence_length = 0;
        (*sequences)[suffixes[i].sid].get_sequence(sequence_data, sequence_length);
        if (suffixes[i].pid + suffixes[i].k <= sequence_length) {
            suffixes[i].bid = bucket_index(
                    &sequence_data[suffixes[i].pid], suffixes[i].k+1);
            suffixes[i].k += 1;
            suffixes[i].next = NULL;
        }
        else {
            /* adding 1 to kmer stepped past the end of the sequence */
            /* invalidate the suffix */
            suffixes[i].sid = npos;
            suffixes[i].pid = npos;
            suffixes[i].bid = npos;
            suffixes[i].k = -1;
            suffixes[i].next = NULL;
        }
    }
    
    /* re-sort the suffixes now that the bucket IDs are updated */
    ::std::sort(local_suffixes, local_suffixes+local_suffixes_size,
            SuffixBucketCompare);

    /* reindex suffixes and add new buckets to end */
    for (size_t i=0; i<local_suffixes_size; ++i) {
        /* skip invalid suffixes */
        if (suffixes[i].bid == npos) {
            continue;
        }
        /* first suffix or current suffix doesn't match previous */
        if (i == 0 || suffixes[i].bid != suffixes[i-1].bid) {
            buckets.push_back(BucketMeta(
                        bucket.offset + i, 1, suffixes[i].bid, suffixes[i].k));
        }
        else if (suffixes[i].bid == suffixes[i-1].bid) {
            suffixes[i-1].next = &suffixes[i];
            buckets.back().size += 1;
        }
        else {
            assert(0); /* this shouldn't be reached */
        }
    }

    /* invalidate bucket */
    bucket.offset = npos;
    bucket.size = 0;
    bucket.bid = npos;
    bucket.k = -1;
}

}; /* namespace pgraph */

