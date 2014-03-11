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


PairCheckGlobal::PairCheckGlobal()
    :   PairCheck()
    ,   mutex()
    ,   alloc_id()
    ,   am_handle()
    ,   s_pairs()
    ,   v_pairs()
{
    am_handle = theAm().amRegister(
            amLocalClient,
            amPostPutServer,
            NULL/*amRemoteClient*/,
            NULL/*amLocalServer*/,
            amPrePutServer);
    alloc_id = theRma().allocColl(sizeof(SetPair*));
    SetPair **ptr = reinterpret_cast<SetPair**>(
            theRma().lookupPointer(RmaPtr(alloc_id)));
    *ptr = &s_pairs;
    theRegistry().addEntry(NEXT_ID, this);
}


PairCheckGlobal::~PairCheckGlobal()
{
}


bool PairCheckGlobal::send_check_message(const pair<size_t,size_t> &_pair)
{
    cout << "PairCheckGlobal::send_check_message("
        << _pair.first
        << ","
        << _pair.second
        << ")" << endl;
    bool payload = false;
    RmaRequest *localReq = RmaRequest::construct();
    RmaRequest *remoteReq = RmaRequest::construct();
    size_t spair[2] = {_pair.first,_pair.second};
    unsigned long kpair[2] = {_pair.first,_pair.second};
    unsigned long index = k_combination2_inv(kpair);
    ProcRank owner = index % (theTwoSided().numProcRanks().toInt());
    PairCheckArg &arg = *(new PairCheckArg(RmaPtr(alloc_id)));
    arg.data = spair;
    arg.dlen = 2*sizeof(size_t);
    arg.rdata = &payload;
    arg.rdlen = sizeof(bool);
    theAm().amGet(owner, arg, am_handle, localReq, remoteReq);
    while (!localReq->test() || !remoteReq->test()) {
        AmListenObjCodelet<NullMutex>* lcodelet;
        if ((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
            lcodelet->execute();
        }
    }
    return payload;
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


void PairCheckGlobal::amLocalClient(const AmContext * const context)
{
    cout << "PairCheckGlobal::amLocalClient" << endl;
    PairCheckArg& arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    delete &arg;
}


bool PairCheckGlobal::do_check(size_t s_pair[2])
{
    LockGuard<PthreadMutex> guard(mutex);
    pair<SetPair::iterator,bool> result;
    result = s_pairs.insert(make_pair(s_pair[0],s_pair[1]));
    return result.second;
}


void PairCheckGlobal::amPostPutServer(const AmContext * const context)
{
    PairCheckArg& arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    SetPair **ptr = reinterpret_cast<SetPair**>(
            theRma().lookupPointer(arg.ptr));
    SetPair *the_set = *ptr;
    size_t *kpair = reinterpret_cast<size_t*>(arg.data);
    PairCheckGlobal *checker = theRegistry().lookup<PairCheckGlobal*>(NEXT_ID);
    bool result = checker->do_check(kpair);
    cout << "PairCheckGlobal::amPostPutServer "
        << kpair[0] << "," << kpair[1]
        << "=" << result << endl;
    theAm().amResponse(*context->arg, &result, sizeof(bool));
    // Data points to a temporary buffer
    delete [] reinterpret_cast<char*>(arg.data);
}


void PairCheckGlobal::amRemoteClient(const AmContext * const context)
{
    cout << "PairCheckGlobal::amRemoteClient" << endl;
}


void PairCheckGlobal::amLocalServer(const AmContext * const context)
{
    cout << "PairCheckGlobal::amLocalServer" << endl;
}


void PairCheckGlobal::amPrePutServer(const AmContext * const context)
{
    cout << "PairCheckGlobal::amPrePutServer" << endl;
    PairCheckArg& arg = *reinterpret_cast<PairCheckArg*>(context->arg);
    // Store into temporary buffer to apply operator later
    assert(arg.dlen.toInt() == 2*sizeof(unsigned long));
    theAm().amResponsePrePut(*context->arg,
            new char[arg.dlen.toInt()], arg.dlen);
}

}; /* namespace pgraph */
