#include <utility>

#include <tascel.h>

#include <iostream>
using ::std::cout;
using ::std::endl;

#include "combinations.h"
#include "PairCheckGlobal.hpp"

using ::std::make_pair;
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

    PairCheckArg()
        :   queueID(-1)
        ,   proc(-1)
        ,   thd(-1)
        ,   id1(-1)
        ,   id2(-1)
        ,   answer(false)
    {}
};



PairCheckGlobal::PairCheckGlobal(int thd)
    :   PairCheck()
    ,   mutex()
    ,   thd(thd)
    ,   alloc_id()
    ,   try_check()
    ,   check_complete()
    ,   s_pairs()
    ,   v_pairs()
    ,   check_answer(false)
    ,   check_response(false)
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
    theRegistry().addEntry(NEXT_ID + thd, this);
}


PairCheckGlobal::~PairCheckGlobal()
{
}


bool PairCheckGlobal::send_check_message(const pair<size_t,size_t> &_pair)
{
#if 0
    cout << "PairCheckGlobal::send_check_message("
        << _pair.first
        << ","
        << _pair.second
        << ")" << endl;
#endif

    assert(check_response == false);

    PairCheckArg msg;
    {
#if defined(THREADED)
        LockGuard<PThreadMutex> guard(mutex);
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
    check_response = false;

    {
#if defined(THREADED)
        LockGuard<PThreadMutex> guard(mutex);
#endif
        return check_answer;
    }
}


PairCheck::SetPair PairCheckGlobal::check(const PairCheck::SetPair &new_pairs)
{
    PairCheck::SetPair ret;

    for (PairCheck::SetPair::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        if (send_check_message(*it)) {
            ret.insert(*it);
        }
    }
    
    return ret;
}


PairCheck::VecPair PairCheckGlobal::check(const PairCheck::VecPair &new_pairs)
{
    assert(0);
    return PairCheck::VecPair();
}


void PairCheckGlobal::try_check_function(const AmContext * const context)
{
    PairCheckArg &arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    PairCheckGlobal &checker = *theRegistry().lookup<PairCheckGlobal*>(arg.queueID);

    {
#if defined(THREADED)
        LockGuard<PThreadMutex> guard(checker.mutex);
#endif
        PairCheckArg &msg = *new PairCheckArg;
        AmRequest *lReq = AmRequest::construct();
        msg.queueID = arg.queueID;
        msg.proc = arg.proc;
        msg.thd = arg.thd;
        msg.id1 = arg.id1;
        msg.id2 = arg.id2;
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
    checker.check_answer = arg.answer;
    checker.check_response = true;
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

}; /* namespace pgraph */
