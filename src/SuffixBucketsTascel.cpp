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

#include "mpix.hpp"
#include "mpix_types.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixBucketsTascel.hpp"

using namespace tascel;
using ::std::cout;
using ::std::endl;
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
    double time = 0.0;

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

    time = MPI_Wtime();
    size_t initial_suffixes_size = 0;
    for (size_t i = start; i < stop; ++i) {
        size_t sequence_length = sequences->get_sequence_size(i);
        if (sequence_length >= size_t(param.window_size)) {
            initial_suffixes_size += sequence_length - param.window_size + 1;
        }
    }
    time = MPI_Wtime() - time;
    vector<Suffix> initial_suffixes(initial_suffixes_size);

#if DEBUG || 1
    {
        vector<double> times = mpix::gather(time, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<times.size(); ++i) {
                stats.push_back(times[i]);
            }
            cout << "                 " << Stats::header() << endl;
            cout << "suffix size time " << stats << endl;
        }
    }
#endif

    /* slide k-mers for every sequence and bucket them */

#if DEBUG || 1
    mpix::print_sync("initial_suffixes_size", initial_suffixes_size, comm);

    {
        vector<size_t> sizes = mpix::gather(initial_suffixes_size, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<sizes.size(); ++i) {
                stats.push_back(sizes[i]);
            }
            cout << "                      " << Stats::header() << endl;
            cout << "initial_suffixes_size " << stats << endl;
        }
    }
#endif

    /* slide k-mers for every sequence and bucket them */
    time = MPI_Wtime();
    for (size_t i = start; i < stop; ++i) {
        size_t stop_index = 0;
        const char *sequence_data = NULL;
        size_t sequence_length = 0;
        Sequence *sequence = NULL;

        sequence = sequences->get_sequence(i);
        sequence->get_sequence(sequence_data, sequence_length);
        if (sequence_length <= ((unsigned)param.window_size)) {
            n_suffixes_skipped += 1;
            continue;
        }
        /* stop_index stops before the assumed DOLLAR terminal character */
        stop_index = sequence_length - param.window_size - 1;
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
        delete sequence;
    }
    initial_suffixes.resize(suffix_index);
    time = MPI_Wtime() - time;
    {
        vector<double> times = mpix::gather(time, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<times.size(); ++i) {
                stats.push_back(times[i]);
            }
            cout << "           " << Stats::header() << endl;
            cout << "k-mer time " << stats << endl;
        }
    }
#if DEBUG || 1
    mpix::print_sync("n_suffixes_skipped", n_suffixes_skipped, comm);
    {
        vector<size_t> skipped = mpix::gather(n_suffixes_skipped, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<skipped.size(); ++i) {
                stats.push_back(skipped[i]);
            }
            cout << "                   " << Stats::header() << endl;
            cout << "n_suffixes_skipped " << stats << endl;
        }
    }
#endif

#if DEBUG || 1
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

#if DEBUG || 1
    mpix::print_zero("before amount_to_send all to all", comm);
#endif
    mpix::alltoall(amount_to_send, amount_to_recv, comm);
#if DEBUG || 1
    mpix::print_zero("after amount_to_send all to all", comm);
#endif

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

#if DEBUG || 1
    mpix::print_zero("before suffix sort", comm);
#endif
    /* We are preparing for the all to all, so we sort the suffixes.
     * Sort is based on owner ID. */
    time = MPI_Wtime();
    std::sort(initial_suffixes.begin(),
              initial_suffixes.end(),
              SuffixOwnerCompareFunctor(comm_size)
    );
    time = MPI_Wtime() - time;
    {
        vector<double> times = mpix::gather(time, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<times.size(); ++i) {
                stats.push_back(times[i]);
            }
            cout << "                 " << Stats::header() << endl;
            cout << "suffix sort time " << stats << endl;
        }
    }
#if DEBUG || 1
    mpix::print_zero("after suffix sort", comm);
#endif

    aid_suffixes = theRma().allocColl(sizeof(Suffix)*suffixes_size);
    suffixes = reinterpret_cast<Suffix*>(
            theRma().lookupPointer(RmaPtr(aid_suffixes)));
#if DEBUG || 1
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
    MPI_Datatype SuffixType = mpix::get_mpi_datatype(initial_suffixes[0]);
#if DEBUG
    mpix::print_sync("initial_suffixes", initial_suffixes, comm);
#endif

#if DEBUG || 1
    mpix::print_zero("before suffix all to all", comm);
#endif
    ierr = MPI_Alltoallv(
            &initial_suffixes[0], &amount_to_send[0],
            &send_displacements[0], SuffixType,
            &suffixes[0], &amount_to_recv[0],
            &recv_displacements[0], SuffixType, comm);
    assert(MPI_SUCCESS == ierr);
#if DEBUG || 1
    mpix::print_zero("after suffix all to all", comm);
#endif

#if DEBUG
    mpix::print_sync("suffixes", suffixes, suffixes_size, comm);
#endif

    vector<BucketMeta> initial_buckets;

    time = MPI_Wtime();
    if (suffixes_size > 0) {
        /* sort the received Suffix instances based on bucket ID */
        std::sort(suffixes, suffixes+suffixes_size, SuffixBucketCompare);
        /* The suffixes contains sorted suffixes based on the buckets
         * they belong to. That means we can simply update the 'next'
         * links! */
        size_t last_id = sequences->size();
        size_t longest = sequences->longest();
        for (size_t i=0; i<suffixes_size; ++i) {
            assert(suffixes[i].sid < last_id);
            assert(suffixes[i].pid < longest);
            assert(suffixes[i].bid < n_buckets);
            assert(suffixes[i].bid % comm_size == comm_rank);
            assert(suffixes[i].k == param.window_size);
            assert(suffixes[i].next == NULL);
            if (i == 0 || suffixes[i-1].bid != suffixes[i].bid) {
                initial_buckets.push_back(BucketMeta(
                            i, 1, suffixes[i].bid, suffixes[i].k));
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
    time = MPI_Wtime() - time;
    {
        vector<double> times = mpix::gather(time, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<times.size(); ++i) {
                stats.push_back(times[i]);
            }
            cout << "                         " << Stats::header() << endl;
            cout << "BucketMeta creation time " << stats << endl;
        }
    }

    /* now that initial bucket are created, collect sizes stats */
    for (size_t i=0; i<initial_buckets.size(); ++i) {
        buckets_stats.push_back(initial_buckets[i].size);
#if DEBUG
        cout << "initial_buckets[" << i << "]=(" << initial_buckets[i].offset
            << "," << initial_buckets[i].size
            << "," << initial_buckets[i].bid
            << "," << initial_buckets[i].k
            << ")" << endl;
#endif
    }
    {
        Stats *all_stats = new Stats[comm_size];
        mpix::gather(&buckets_stats, 1, all_stats, 1, 0, comm);
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
        if (0 == comm_rank) {
            cout << buckets_stats << endl;
        }
        delete [] all_stats;
    }

#if DEBUG || 1
    {
        Stats *all_stats = new Stats[comm_size];
        mpix::gather(&buckets_stats, 1, all_stats, 1, 0, comm);
        if (0 == comm_rank) {
            cout << Stats::header() << endl;
            for (size_t i=0; i<comm_size; ++i) {
                cout << all_stats[i] << endl;
            }
        }
        delete [] all_stats;
    }
#endif

#if DEBUG || 1
    mpix::print_sync("buckets_size_before", initial_buckets.size(), comm);
#endif

    /* now that initial buckets are created, refine the big ones */
    time = MPI_Wtime();
    bool done = false;
    start = 0;
    while (!done) {
        done = true;
        vector<BucketMeta> new_buckets;
        for (size_t i=start; i<initial_buckets.size(); ++i) {
            BucketMeta &bucket = initial_buckets[i];
#if DEBUG
            if (0 == comm_rank) {
                ::std::cout << "bucket.size >? bucket_stats.mean() + bucket_cutoff*buckets_stats.stddev() -- "
                    << bucket.size << " >? " << buckets_stats.mean() + param.bucket_cutoff*buckets_stats.stddev()
                    << ::std::endl;
                ::std::cout << "bucket.k <? param.exact_match_length -- "
                    << bucket.k << " <? " << param.exact_match_length
                    << ::std::endl;
            }
#endif
            if (bucket.size > buckets_stats.mean() + param.bucket_cutoff*buckets_stats.stddev()
                    && bucket.k+1 < param.exact_match_length) {
                vector<BucketMeta> local_new_buckets;
                refine_bucket(bucket, local_new_buckets);
                done = false;
#if DEBUG
                ::std::cout << "created " << local_new_buckets.size()
                    << " new buckets" << ::std::endl;
#endif
                new_buckets.insert(new_buckets.end(),
                        local_new_buckets.begin(), local_new_buckets.end());
                /* invalidate bucket since it was refined */
                bucket.offset = npos;
                bucket.size = 0;
                bucket.bid = npos;
                bucket.k = -1;
            }
        }
        start = initial_buckets.size();
        initial_buckets.insert(initial_buckets.end(),
                new_buckets.begin(), new_buckets.end());
    }
    buckets_size = initial_buckets.size();
    time = MPI_Wtime() - time;
    {
        vector<double> times = mpix::gather(time, 0, comm);
        if (0 == comm_rank) {
            Stats stats;
            for (size_t i=0; i<times.size(); ++i) {
                stats.push_back(times[i]);
            }
            cout << "                " << Stats::header() << endl;
            cout << "refinement time " << stats << endl;
        }
    }

#if DEBUG || 1
    mpix::print_sync("buckets_size", buckets_size, comm);
#endif
    aid_meta = theRma().allocColl(sizeof(BucketMeta)*buckets_size);
    buckets = reinterpret_cast<BucketMeta*>(
            theRma().lookupPointer(RmaPtr(aid_meta)));
    ::std::copy(initial_buckets.begin(), initial_buckets.end(), buckets);
#if DEBUG
    for (size_t i=0; i<initial_buckets.size(); ++i) {
        buckets[i] = initial_buckets[i];
        cout << "initial_buckets[" << i << "]=(" << initial_buckets[i].offset
            << "," << initial_buckets[i].size
            << "," << initial_buckets[i].bid
            << "," << initial_buckets[i].k
            << ")" << endl;
        cout << "buckets[" << i << "]=(" << buckets[i].offset
            << "," << buckets[i].size
            << "," << buckets[i].bid
            << "," << buckets[i].k
            << ")" << endl;
    }
#endif

    buckets_size_global = buckets_size;
    mpix::allreduce(buckets_size_global, MPI_SUM, comm);
#if DEBUG || 1
    mpix::print_zero("buckets_size_global", buckets_size_global, comm);
#endif

    /* now that bucket are refined, re-collect sizes stats */
    buckets_stats.clear();
    for (size_t i=0; i<initial_buckets.size(); ++i) {
        buckets_stats.push_back(initial_buckets[i].size);
    }
    {
        Stats *all_stats = new Stats[comm_size];
        mpix::gather(&buckets_stats, 1, all_stats, 1, 0, comm);
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
        if (0 == comm_rank) {
            cout << buckets_stats << endl;
        }
        delete [] all_stats;
    }

#if DEBUG || 1
    {
        Stats *all_stats = new Stats[comm_size];
        mpix::gather(&buckets_stats, 1, all_stats, 1, 0, comm);
        if (0 == comm_rank) {
            cout << Stats::header() << endl;
            for (size_t i=0; i<comm_size; ++i) {
                cout << all_stats[i] << endl;
            }
        }
        delete [] all_stats;
    }
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

    //cout << "getting bucket (" << owner << "," << bucket_index << ")" << endl;
    if (owner == int(comm_rank)) {
        /* already owned, just return it */
        BucketMeta &meta = buckets[bucket_index];
#if DEBUG
        cout << "owned bucket = (" << meta.offset
            << "," << meta.size
            << "," << meta.bid
            << "," << meta.k << ")" << endl;
#endif
        bucket->size = meta.size;
        bucket->bid = meta.bid;
        bucket->k = meta.k;
        bucket->owner = owner;
        //cout << "owner, size was " << bucket->size << endl;
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
    if (bucket->owner != int(comm_rank)
            && bucket->size > 0
            && bucket->suffixes != NULL) {
        delete [] bucket->suffixes;
    }
    delete bucket;
}


void SuffixBucketsTascel::refine_bucket(
        BucketMeta bucket,
        vector<BucketMeta> &result)
{
    Suffix *local_suffixes = &suffixes[bucket.offset];
    size_t local_suffixes_size = bucket.size;

#if DEBUG
    ::std::cout << "refining bucket"
        << " " << bucket.bid
        << " " << bucket.k
        << " " << bucket_kmer(bucket.bid, bucket.k)
        << " " << bucket.offset
        << " " << bucket.size
        << ::std::endl;
#endif

    for (size_t i=0; i<local_suffixes_size; ++i) {
        const char *sequence_data = NULL;
        size_t sequence_length = 0;
        Sequence *sequence = NULL;
        //cout << "suffixes[" << bucket.offset+i << "]=" << local_suffixes[i] << endl;
        if (local_suffixes[i].bid != bucket.bid) {
            cout << "suffixes[" << bucket.offset+i << "].bid=" << local_suffixes[i].bid
                << " != bucket.bid=" << bucket.bid << endl;
            MPI_Abort(comm, -1);
        }
        assert(local_suffixes[i].bid == bucket.bid);
        assert(local_suffixes[i].bid != npos);
        sequence = sequences->get_sequence(local_suffixes[i].sid);
        sequence->get_sequence(sequence_data, sequence_length);
        if (local_suffixes[i].pid + local_suffixes[i].k + 1 < sequence_length) {
            size_t new_bid = bucket_index(
                    &sequence_data[local_suffixes[i].pid], local_suffixes[i].k+1);
#if DEBUG
            cout << "suffixes[" << bucket.offset+i << "].bid = " << local_suffixes[i].bid
                << " new_bid = " << new_bid
                << endl;
#endif
            local_suffixes[i].bid = new_bid;
            local_suffixes[i].k += 1;
            local_suffixes[i].next = NULL;
        }
        else {
#if DEBUG
            cout << "suffixes[" << bucket.offset+i << "].bid = " << local_suffixes[i].bid
                << " stepped off end of sequence" << endl;
#endif
            /* adding 1 to kmer stepped past the end of the sequence */
            /* invalidate the suffix */
            local_suffixes[i].sid = npos;
            local_suffixes[i].pid = npos;
            local_suffixes[i].bid = npos;
            local_suffixes[i].k = -1;
            local_suffixes[i].next = NULL;
        }
        delete sequence;
    }
    
    /* re-sort the suffixes now that the bucket IDs are updated */
    ::std::sort(local_suffixes, local_suffixes+local_suffixes_size,
            SuffixBucketCompare);

    /* reindex suffixes and add new buckets to end */
    for (size_t i=0; i<local_suffixes_size; ++i) {
        /* skip invalid suffixes */
        if (local_suffixes[i].bid == npos) {
            continue;
        }
        /* first suffix or current suffix doesn't match previous */
        if (i == 0 || local_suffixes[i].bid != local_suffixes[i-1].bid) {
            size_t offset_new = bucket.offset + i;
            size_t size_new = 1;
            size_t bid_new = local_suffixes[i].bid;
            int k_new = local_suffixes[i].k;
            result.push_back(BucketMeta(offset_new, size_new, bid_new, k_new));
                        //bucket.offset + i, 1, local_suffixes[i].bid, local_suffixes[i].k));
        }
        else if (local_suffixes[i].bid == local_suffixes[i-1].bid) {
            local_suffixes[i-1].next = &local_suffixes[i];
            result.back().size += 1;
        }
        else {
            assert(0); /* this shouldn't be reached */
        }
    }
}

}; /* namespace pgraph */

