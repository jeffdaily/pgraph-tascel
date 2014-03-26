/**
 * @file PairCheckGlobalServer.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <tascel.h>

//#include <algorithm>
#include <map>

#include "combinations.h"
#include "PairCheckGlobalServer.hpp"

using ::std::map;
using namespace tascel;

namespace pgraph {


#define NEXT_ID (747274)


#define DEBUG 0
#if DEBUG
#include <iostream>
using ::std::cout;
using ::std::endl;
static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}
#endif


struct BulkPairCheckArg : public AmArg {
    int queueID;
    int proc;
    int thd;
    size_t last_size;
#if USE_SET
    SetPair *pairs;
#else
    VecPair *pairs;
#endif

    BulkPairCheckArg()
        :   queueID(-1)
        ,   proc(-1)
        ,   thd(-1)
        ,   pairs(NULL)
    {}
};


PairCheckGlobalServer::PairCheckGlobalServer(int thd)
    :   PairCheck()
    ,   mutex()
    ,   thd(thd)
    ,   bulk_try_check()
    ,   bulk_check_complete()
    ,   dispatcher()
    ,   last_size(0)
    ,   server_response(0)
{
    bulk_try_check = theAm().amRegister(
            bulk_try_check_local_function,
            bulk_try_check_function,
            NULL,
            NULL,
            NULL);
    bulk_check_complete = theAm().amRegister(
            bulk_check_complete_local_function,
            bulk_check_complete_function,
            NULL,
            NULL,
            NULL);
    theRegistry().addEntry(NEXT_ID + thd, this);
}


PairCheckGlobalServer::~PairCheckGlobalServer()
{
}


#if USE_SET
SetPair PairCheckGlobalServer::check(const SetPair &new_pairs)
{
    SetPair ret;
    int proc = theTwoSided().getProcRank().toInt();
    int nprocs = theTwoSided().numProcRanks().toInt();

#if DEBUG
    cout << "PairCheckGlobalServer::check begin"
        << " rank=" << proc
        << " thd=" << thd
        << " trank=" << trank(thd)
        << endl;
#endif
    assert(server_response == 0);

    /* sort the pairs into buckets based on their worker rank */
    map<int,VecPair> pairs_parted;
    for (SetPair::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        unsigned long k_pair[2] = {it->first,it->second};
        unsigned long index = k_combination2_inv(k_pair);
        int owner = index % nprocs;
        pairs_parted[owner].push_back(*it);
    }

    map<int,VecPair>::iterator it;
    for (it=pairs_parted.begin(); it!=pairs_parted.end(); ++it) {
        BulkPairCheckArg &msg = *new BulkPairCheckArg;
        const int &owner = it->first;
        VecPair &the_pairs = it->second;
        //Pair *copy_pairs = new Pair[the_pairs.size()];
        //::std::copy(the_pairs.begin(), the_pairs.end(), copy_pairs);
        AmRequest *lReq = AmRequest::construct();
        msg.queueID = owner % NUM_WORKERS;
        msg.proc = proc;
        msg.thd = thd;
        //msg.data = copy_pairs;
        msg.data = &the_pairs[0];
        msg.dlen = sizeof(Pair)*the_pairs.size();
        msg.pairs = &ret;
#if defined(THREADED)
        {
            LockGuard<PthreadMutex> guard(mutex);
            server_response += 1;
        }
#else
        server_response += 1;
#endif
        theAm().amPut(owner, msg, bulk_try_check, lReq);
        dispatcher.registerCodelet(lReq);
#if DEBUG
        cout << "PairCheckGlobalServer::check send"
            << " rank=" << proc
            << " thd=" << thd
            << " trank=" << trank(thd)
            << " server_response=" << server_response
            << " dest=" << owner/NUM_WORKERS
            << " queueID=" << msg.queueID
            << endl;
#endif
    }

#if DEBUG
    cout << "PairCheckGlobalServer::check waiting"
        << " rank=" << proc
        << " thd=" << thd
        << " trank=" << trank(thd)
        << endl;
#endif
    while (server_response>0) {
#if !defined(THREADED)
        AmListenObjCodelet<NullMutex>* listenCodelet;
        if ((listenCodelet = theAm().amListeners[0]->progress()) != NULL) {
            listenCodelet->execute();
        }
#endif
    }
    while (!dispatcher.empty()) {
        Codelet* codelet;
        if ((codelet = dispatcher.progress()) != NULL) {
            codelet->execute();
            delete reinterpret_cast<AmRequest*>(codelet);
        }
    }
    assert(server_response == 0);
    assert(dispatcher.empty());
    
#if DEBUG
    cout << "PairCheckGlobalServer::check end"
        << " rank=" << proc
        << " thd=" << thd
        << " trank=" << trank(thd)
        << endl;
#endif
    return ret;
}
#else
SetPair PairCheckGlobalServer::check(const SetPair &new_pairs)
{
    VecPair new_pairs_as_vec(new_pairs.begin(), new_pairs.end());
    VecPair ret_pairs_as_vec = check(new_pairs_as_vec);
    return SetPair(ret_pairs_as_vec.begin(), ret_pairs_as_vec.end());
}
#endif


#if USE_SET
VecPair PairCheckGlobalServer::check(const VecPair &new_pairs)
{
    SetPair new_pairs_as_set(new_pairs.begin(), new_pairs.end());
    SetPair ret_pairs_as_set = check(new_pairs_as_set);
    return VecPair(ret_pairs_as_set.begin(), ret_pairs_as_set.end());
}
#else
VecPair PairCheckGlobalServer::check(const VecPair &new_pairs)
{
    VecPair ret;
    int proc = theTwoSided().getProcRank().toInt();
    int nprocs = theTwoSided().numProcRanks().toInt();

#if DEBUG
    cout << "PairCheckGlobalServer::check begin"
        << " rank=" << proc
        << " thd=" << thd
        << " trank=" << trank(thd)
        << endl;
#endif
    assert(server_response == 0);

    /* sort the pairs into buckets based on their worker rank */
    map<int,VecPair> pairs_parted;
    for (VecPair::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        unsigned long k_pair[2] = {it->first,it->second};
        unsigned long index = k_combination2_inv(k_pair);
        int owner = index % nprocs;
        pairs_parted[owner].push_back(*it);
    }

    map<int,VecPair>::iterator it;
    for (it=pairs_parted.begin(); it!=pairs_parted.end(); ++it) {
        BulkPairCheckArg &msg = *new BulkPairCheckArg;
        const int &owner = it->first;
        VecPair &the_pairs = it->second;
        //Pair *copy_pairs = new Pair[the_pairs.size()];
        //::std::copy(the_pairs.begin(), the_pairs.end(), copy_pairs);
        AmRequest *lReq = AmRequest::construct();
        msg.queueID = owner % NUM_WORKERS;
        msg.proc = proc;
        msg.thd = thd;
        //msg.data = copy_pairs;
        msg.data = &the_pairs[0];
        msg.dlen = sizeof(Pair)*the_pairs.size();
        msg.pairs = &ret;
#if defined(THREADED)
        {
            LockGuard<PthreadMutex> guard(mutex);
            server_response += 1;
        }
#else
        server_response += 1;
#endif
        theAm().amPut(owner, msg, bulk_try_check, lReq);
        dispatcher.registerCodelet(lReq);
#if DEBUG
        cout << "PairCheckGlobalServer::check send"
            << " rank=" << proc
            << " thd=" << thd
            << " trank=" << trank(thd)
            << " server_response=" << server_response
            << " dest=" << owner/NUM_WORKERS
            << " queueID=" << msg.queueID
            << endl;
#endif
    }

#if DEBUG
    cout << "PairCheckGlobalServer::check waiting"
        << " rank=" << proc
        << " thd=" << thd
        << " trank=" << trank(thd)
        << endl;
#endif
    while (server_response>0) {
#if !defined(THREADED)
        AmListenObjCodelet<NullMutex>* listenCodelet;
        if ((listenCodelet = theAm().amListeners[0]->progress()) != NULL) {
            listenCodelet->execute();
        }
#endif
    }
    while (!dispatcher.empty()) {
        Codelet* codelet;
        if ((codelet = dispatcher.progress()) != NULL) {
            codelet->execute();
            delete reinterpret_cast<AmRequest*>(codelet);
        }
    }
    assert(server_response == 0);
    assert(dispatcher.empty());
    
#if DEBUG
    cout << "PairCheckGlobalServer::check end"
        << " rank=" << proc
        << " thd=" << thd
        << " trank=" << trank(thd)
        << endl;
#endif
    return ret;
}
#endif


size_t PairCheckGlobalServer::size()
{
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(mutex);
#endif
    return last_size;
}


void PairCheckGlobalServer::bulk_try_check_function(const AmContext * const context)
{
    static SetPair s_pairs;
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
    PairCheckGlobalServer &checker = *theRegistry().lookup<PairCheckGlobalServer*>(NEXT_ID+arg.queueID);
#if DEBUG
    cout << "PairCheckGlobalServer::bulk_try_check_function"
        << " trank=" << trank(checker.thd)
        << " arg.proc=" << arg.proc
        << " arg.thd=" << arg.thd
        << endl;
#endif

    BulkPairCheckArg &msg = *new BulkPairCheckArg;
    AmRequest *lReq = AmRequest::construct();
    msg.queueID = arg.thd;
    msg.proc = arg.proc;
    msg.thd = arg.thd;
    msg.pairs = arg.pairs;

    Pair *pairs_to_check = reinterpret_cast<Pair*>(arg.data);
    int n = arg.dlen.toInt() / sizeof(Pair);
    assert(arg.dlen % sizeof(Pair) == 0);
    Pair *pairs_to_return = new Pair[n];
    int size = 0;
    for (int i=0; i<n; ++i) {
        pair<SetPair::iterator,bool> result =
            s_pairs.insert(pairs_to_check[i]);
        if (result.second) {
            pairs_to_return[size++] = pairs_to_check[i];
        }
    }
    msg.last_size = s_pairs.size();
    if (size > 0) {
        msg.data = pairs_to_return;
        msg.dlen = sizeof(Pair) * size;
    }
    else {
        delete [] pairs_to_return;
        msg.data = NullPtr;
        msg.dlen = 0;
    }
    theAm().amPut(arg.proc, msg, checker.bulk_check_complete, lReq);
#if defined(THREADED)
    serverDispatcher.registerCodelet(lReq);
#else
    checker.dispatcher.registerCodelet(lReq);
#endif
}


void PairCheckGlobalServer::bulk_try_check_local_function(const AmContext * const context)
{
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
#if DEBUG
    PairCheckGlobalServer &checker = *theRegistry().lookup<PairCheckGlobalServer*>(NEXT_ID+arg.queueID);
    cout << "PairCheckGlobalServer::bulk_try_check_local_function " << trank(checker.thd) << endl;
#endif
    /* don't delete the data since it is internal to a vector */
#if 0
    if (arg.data != NullPtr) {
        assert(arg.dlen > 0);
        Pair *pairs_that_were_sent_to_check = reinterpret_cast<Pair*>(arg.data);
        delete [] pairs_that_were_sent_to_check;
        arg.data = NullPtr;
    }
    else {
        assert(arg.data == NullPtr);
        assert(arg.dlen == 0);
    }
#endif
    delete &arg;
}


void PairCheckGlobalServer::bulk_check_complete_function(const AmContext * const context)
{
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
    PairCheckGlobalServer &checker = *theRegistry().lookup<PairCheckGlobalServer*>(NEXT_ID+arg.queueID);
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(checker.mutex);
#endif
#if DEBUG
    cout << "PairCheckGlobalServer::bulk_check_complete_function " << trank(checker.thd) << endl;
#endif
    assert(checker.server_response >= 1);
    Pair *pairs_to_insert = reinterpret_cast<Pair*>(arg.data);
    int n = arg.dlen.toInt() / sizeof(Pair);
    assert(arg.dlen % sizeof(Pair) == 0);
#if USE_SET
    arg.pairs->insert(pairs_to_insert, pairs_to_insert+n);
#else
    arg.pairs->insert(arg.pairs->end(), pairs_to_insert, pairs_to_insert+n);
#endif
    checker.last_size = arg.last_size;
    checker.server_response -= 1;
#if DEBUG
    cout << "PairCheckGlobalServer::bulk_check_complete_function " << trank(checker.thd)
        << " server_response=" << checker.server_response << endl;
#endif
    assert(checker.server_response >= 0);
}


void PairCheckGlobalServer::bulk_check_complete_local_function(const AmContext * const context)
{
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
#if DEBUG
    PairCheckGlobalServer &checker = *theRegistry().lookup<PairCheckGlobalServer*>(NEXT_ID+arg.queueID);
    cout << "PairCheckGlobalServer::bulk_check_complete_local_function " << trank(checker.thd) << endl;
#endif
    if (arg.data != NullPtr) {
        assert(arg.dlen > 0);
        Pair *pairs_to_return = reinterpret_cast<Pair*>(arg.data);
        delete [] pairs_to_return;
    }
    else {
        assert(arg.data == NullPtr);
        assert(arg.dlen == 0);
    }
    delete &arg;
}

}; /* namespace pgraph */
