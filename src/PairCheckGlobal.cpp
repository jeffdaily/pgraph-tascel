#include <utility>

#include <tascel.h>

#include <iostream>
#include <map>

#include "combinations.h"
#include "PairCheckGlobal.hpp"

using ::std::cout;
using ::std::endl;
using ::std::make_pair;
using ::std::map;
using namespace tascel;

namespace pgraph {


#define NEXT_ID (747274)


struct PairCheckArg : public AmArg {
    int queueID;
    int proc;
    int thd;
    size_t id1;
    size_t id2;
    bool answer;
    bool *ptr_answer;
    bool *ptr_response;

    PairCheckArg()
        :   queueID(-1)
        ,   proc(-1)
        ,   thd(-1)
        ,   id1(-1)
        ,   id2(-1)
        ,   answer(false)
        ,   ptr_answer(NULL)
        ,   ptr_response(NULL)
    {}
};


struct BulkPairCheckArg : public AmArg {
    int queueID;
    int proc;
    int thd;
    SetPair *ptr_pairs;
    volatile int *ptr_response;
    VecPair *ptr_local_pairs;

    BulkPairCheckArg()
        :   queueID(-1)
        ,   proc(-1)
        ,   thd(-1)
        ,   ptr_pairs(NULL)
        ,   ptr_response(NULL)
        ,   ptr_local_pairs(NULL)
    {}
};


PairCheckGlobal::PairCheckGlobal(int thd)
    :   PairCheck()
    ,   mutex()
    ,   thd(thd)
    ,   alloc_id()
    ,   try_check()
    ,   check_complete()
    ,   bulk_try_check()
    ,   bulk_check_complete()
    ,   s_pairs()
    ,   v_pairs()
    ,   dispatcher()
{
    try_check = theAm().amRegister(
            NULL,
            try_check_function,
            NULL,
            NULL,
            NULL);
    check_complete = theAm().amRegister(
            check_complete_local_function,
            check_complete_function,
            NULL,
            NULL,
            NULL);
    bulk_try_check = theAm().amRegister(
            NULL,
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


PairCheckGlobal::~PairCheckGlobal()
{
}


bool PairCheckGlobal::send_check_message(const pair<size_t,size_t> &_pair)
{
    bool check_response = false;
    bool check_answer = false;

    PairCheckArg msg;
    {
#if defined(THREADED)
        LockGuard<PthreadMutex> guard(mutex);
#endif
        unsigned long k_pair[2] = {_pair.first,_pair.second};
        unsigned long index = k_combination2_inv(k_pair);
        int owner = index % (theTwoSided().numProcRanks().toInt()*NUM_WORKERS);
        AmRequest *lReq = AmRequest::construct();
        msg.queueID = NEXT_ID + owner % NUM_WORKERS;
        msg.proc = theTwoSided().getProcRank().toInt();
        msg.thd = thd;
        msg.id1 = _pair.first;
        msg.id2 = _pair.second;
        msg.answer = false;
        msg.ptr_answer = &check_answer;
        msg.ptr_response = &check_response;
        theAm().amPut(owner / NUM_WORKERS, msg, try_check, lReq);
        dispatcher.registerCodelet(lReq);
    }

    while (!check_response || !dispatcher.empty()) {
        Codelet* codelet;
        if ((codelet = dispatcher.progress()) != NULL) {
            codelet->execute();
            delete reinterpret_cast<AmRequest*>(codelet);
        }
#if !defined(THREADED)
        AmListenObjCodelet<NullMutex>* listenCodelet;
        if ((listenCodelet = theAm().amListeners[0]->progress()) != NULL) {
            listenCodelet->execute();
        }
#endif
    }

    {
#if defined(THREADED)
        LockGuard<PthreadMutex> guard(mutex);
#endif
        return check_answer;
    }
}


void PairCheckGlobal::bulk_send_check_message(
         const SetPair &new_pairs,
         SetPair &result)
{
    volatile int bulk_check_response = 0;

    /* sort the pairs into buckets based on their worker rank */
    int proc = theTwoSided().getProcRank().toInt();
    int nworkers = theTwoSided().numProcRanks().toInt()*NUM_WORKERS;
    map<int,VecPair> pairs_parted;
    for (SetPair::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        unsigned long k_pair[2] = {it->first,it->second};
        unsigned long index = k_combination2_inv(k_pair);
        int owner = index % nworkers;
        pairs_parted[owner].push_back(*it);
    }

    {
#if defined(THREADED)
        LockGuard<PthreadMutex> guard(mutex);
#endif
        map<int,VecPair>::iterator it;
        for (it=pairs_parted.begin(); it!=pairs_parted.end(); ++it) {
            BulkPairCheckArg msg;
            const int &owner = it->first;
            VecPair &the_pairs = it->second;
            AmRequest *lReq = AmRequest::construct();
            msg.queueID = NEXT_ID + owner % NUM_WORKERS;
            msg.proc = proc;
            msg.thd = thd;
            msg.data = &the_pairs[0];
            msg.dlen = sizeof(Pair)*the_pairs.size();
            msg.ptr_pairs = &result;
            msg.ptr_response = &bulk_check_response;
            msg.ptr_local_pairs = NULL;
            bulk_check_response += 1;
            theAm().amPut(owner / NUM_WORKERS, msg, bulk_try_check, lReq);
            dispatcher.registerCodelet(lReq);
        }
    }

    while (bulk_check_response>0 || !dispatcher.empty()) {
        Codelet* codelet;
        if ((codelet = dispatcher.progress()) != NULL) {
            codelet->execute();
            delete reinterpret_cast<AmRequest*>(codelet);
        }
#if !defined(THREADED)
        AmListenObjCodelet<NullMutex>* listenCodelet;
        if ((listenCodelet = theAm().amListeners[0]->progress()) != NULL) {
            listenCodelet->execute();
        }
#endif
    }
    assert(bulk_check_response == 0);
}


SetPair PairCheckGlobal::check(const SetPair &new_pairs)
{
    SetPair ret;

#if 0
    for (SetPair::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        if (send_check_message(*it)) {
            ret.insert(*it);
        }
    }
#else
    bulk_send_check_message(new_pairs, ret);
#endif
    
    return ret;
}


VecPair PairCheckGlobal::check(const VecPair & /*new_pairs*/)
{
    assert(0);
    return VecPair();
}


void PairCheckGlobal::try_check_function(const AmContext * const context)
{
    PairCheckArg &arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);

    {
#if defined(THREADED)
        LockGuard<PthreadMutex> guard(checker.mutex);
#endif
        PairCheckArg &msg = *new PairCheckArg;
        AmRequest *lReq = AmRequest::construct();
        msg.queueID = arg.queueID;
        msg.proc = arg.proc;
        msg.thd = arg.thd;
        msg.id1 = arg.id1;
        msg.id2 = arg.id2;
        msg.ptr_answer = arg.ptr_answer;
        msg.ptr_response = arg.ptr_response;
        pair<SetPair::iterator,bool> result = checker.s_pairs.insert(
                make_pair(arg.id1,arg.id2));
        msg.answer = result.second;
        theAm().amPut(arg.proc, msg, checker.check_complete, lReq);
#if defined(THREADED)
        serverDispatcher.registerCodelet(lReq);
#else
        checker.dispatcher.registerCodelet(lReq);
#endif
    }
}


void PairCheckGlobal::check_complete_function(const AmContext * const context)
{
    PairCheckArg &arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(checker.mutex);
#endif
    *(arg.ptr_answer) = arg.answer;
    *(arg.ptr_response) = true;
}


void PairCheckGlobal::check_complete_local_function(const AmContext * const context)
{
    PairCheckArg &arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(checker.mutex);
#endif
    delete &arg;
}


void PairCheckGlobal::bulk_try_check_function(const AmContext * const context)
{
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);

    {
#if defined(THREADED)
        LockGuard<PthreadMutex> guard(checker.mutex);
#endif
        BulkPairCheckArg &msg = *new BulkPairCheckArg;
        AmRequest *lReq = AmRequest::construct();
        msg.queueID = arg.queueID;
        msg.proc = arg.proc;
        msg.thd = arg.thd;
        msg.ptr_pairs = arg.ptr_pairs;
        msg.ptr_response = arg.ptr_response;
        msg.ptr_local_pairs = new VecPair;

        Pair *pairs_to_check = reinterpret_cast<Pair*>(arg.data);
        int n = arg.dlen.toInt() / sizeof(Pair);
        VecPair &pairs_to_return = *msg.ptr_local_pairs;
        assert(arg.dlen % sizeof(Pair) == 0);
        for (int i=0; i<n; ++i) {
            pair<SetPair::iterator,bool> result =
                checker.s_pairs.insert(pairs_to_check[i]);
            if (result.second) {
                pairs_to_return.push_back(pairs_to_check[i]);
            }
        }
        if (pairs_to_return.size() > 0) {
            msg.data = &pairs_to_return[0];
            msg.dlen = sizeof(Pair) * pairs_to_return.size();
        }
        else {
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
}


void PairCheckGlobal::bulk_check_complete_function(const AmContext * const context)
{
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(checker.mutex);
#endif
    Pair *pairs_to_insert = reinterpret_cast<Pair*>(arg.data);
    int n = arg.dlen.toInt() / sizeof(Pair);
    assert(arg.dlen % sizeof(Pair) == 0);
    (*(arg.ptr_pairs)).insert(pairs_to_insert, pairs_to_insert+n);
    *(arg.ptr_response) -= 1;
}


void PairCheckGlobal::bulk_check_complete_local_function(const AmContext * const context)
{
    BulkPairCheckArg &arg = *reinterpret_cast<BulkPairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(checker.mutex);
#endif
    if (arg.ptr_local_pairs != NullPtr) {
        delete arg.ptr_local_pairs;
    }
    else {
        assert(arg.data == NullPtr);
        assert(arg.dlen == 0);
    }
    delete &arg;
}

}; /* namespace pgraph */
