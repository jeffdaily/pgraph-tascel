#include <utility>

#include <tascel.h>

#include "combinations.h"
#include "PairCheckGlobal.hpp"

using ::std::make_pair;
using namespace tascel;

namespace pgraph {


#define NEXT_ID (747274)


PairCheckGlobal::PairCheckGlobal(int thread_rank)
    :   pairs()
    ,   thread_rank(thread_rank)
    ,   check_id(NEXT_ID + thread_rank)
    ,   try_check()
    ,   complete_check()
    ,   server_mutex()
    ,   dispatcher()
    ,   check_response(0)
{
    theRegistry().addEntry(check_id, this);
    try_check = theAm().amRegister(0, am_try_check_function, 0, 0, 0);
    complete_check = theAm().amRegister(
            am_complete_check_local_function,
            am_complete_check_function, 0, 0, 0);
}


PairCheckGlobal::~PairCheckGlobal()
{
}


void PairCheckGlobal::send_check_message(const pair<size_t,size_t> &_pair)
{
    TryCheck msg;
    AmRequest *lReq = AmRequest::construct();
    unsigned long kpair[] = {_pair.first,_pair.second};
    unsigned long index = k_combination2_inv(kpair);
    ProcRank owner = index % (theTwoSided().numProcRanks().toInt()*NUM_WORKERS);
    msg.id = NEXT_ID + index % NUM_WORKERS;
    msg.proc = theTwoSided().getProcRank().toInt(); /* my rank */
    msg.thd = thread_rank; /* my thread rank */
    msg.a = _pair.first;
    msg.b = _pair.second;
    msg.result = false; /* not used */
    theAm().amPut(owner / NUM_WORKERS, msg, try_check, lReq);
    dispatcher.registerCodelet(lReq);

    {
#if defined(THREADED)
        LockGuard<PthreadMutex> guard(server_mutex);
#endif
        check_response += 1;
    }

    Codelet* codelet;
    if ((codelet = dispatcher.progress()) != NULL) {
        codelet->execute();
        delete reinterpret_cast<AmRequest*>(codelet);
    }
#if !defined(THREADED)
    AmListenObjCodelet<NullMutex> *lcodelet;
    if ((lcodelet = theAm().amListeners[0]->progress()) != NULL) {
        lcodelet->execute();
    }
#endif
}


PairCheck::SetPair PairCheckGlobal::check(const PairCheck::SetPair &new_pairs)
{
    PairCheck::SetPair ret;

    assert(check_response == 0);
    assert(s_pairs_response.empty());
    assert(v_pairs_response.empty());

    for (PairCheck::SetPair::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        send_check_message(*it);
    }

    while (check_response > 0) {
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
    
    return ret;
}


PairCheck::VecPair PairCheckGlobal::check(const PairCheck::VecPair &new_pairs)
{
    assert(check_response == 0);

    return PairCheck::VecPair();
}


void PairCheckGlobal::am_try_check_function(const AmContext * const context)
{
    TryCheck &tc = *reinterpret_cast<TryCheck*>(context->arg);
    PairCheckGlobal &check = *theRegistry().lookup<PairCheckGlobal*>(tc.id);
    check.do_try_check_function(tc);
}


void PairCheckGlobal::do_try_check_function(const TryCheck &tc)
{
#if defined(THREADED)
    LockGuard<PthreadMutex> guard(server_mutex);
#endif

    TryCheck &msg = *new TryCheck;
    AmRequest *lReq = AmRequest::construct();
    msg.proc = tc.proc;
    msg.id = NEXT_ID + tc.thd;
    msg.thd = tc.thd;
    msg.a = tc.a;
    msg.b = tc.b;
    msg.result = pairs.find(make_pair(tc.a,tc.b)) == pairs.end();
    theAm().amPut(tc.proc, msg, complete_check, lReq);
#if defined(THREADED)
    serverDispatcher.registerCodelet(lReq);
#else
    dispatcher.registerCodelet(lReq);
#endif
}


void PairCheckGlobal::am_complete_check_function(const AmContext * const context)
{
    TryCheck &tc = *reinterpret_cast<TryCheck*>(context->arg);
    PairCheckGlobal &check = *theRegistry().lookup<PairCheckGlobal*>(tc.id);
    check.do_complete_check_function(tc);
}


void PairCheckGlobal::do_complete_check_function(const TryCheck &tc)
{
}


void PairCheckGlobal::am_complete_check_local_function(const AmContext * const context)
{
    TryCheck &tc = *reinterpret_cast<TryCheck*>(context->arg);
    PairCheckGlobal &check = *theRegistry().lookup<PairCheckGlobal*>(tc.id);
    check.do_complete_check_local_function(tc);
}


void PairCheckGlobal::do_complete_check_local_function(const TryCheck &tc)
{
}


}; /* namespace pgraph */
