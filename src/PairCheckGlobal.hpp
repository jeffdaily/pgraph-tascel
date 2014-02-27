/**
 * @file PairCheckGlobal.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECKGLOBAL_H_
#define _PGRAPH_PAIRCHECKGLOBAL_H_

#include <tascel.h>

#include "PairCheck.hpp"

using ::tascel::AmArg;
using ::tascel::AmContext;
using ::tascel::AmHandle;
using ::tascel::Dispatcher;
using ::tascel::PthreadMutex;

namespace pgraph {

class PairCheckGlobal : public PairCheck
{
    public:
        PairCheckGlobal(int thread_rank);
        virtual ~PairCheckGlobal();

        virtual SetPair check(const SetPair &pairs);
        virtual VecPair check(const VecPair &pairs);

    protected:
        struct TryCheck : public AmArg {
            int id;
            int proc;
            int thd;
            size_t a;
            size_t b;
            bool result;

            TryCheck()
                :   id(-1)
                    ,   proc(-1)
                    ,   thd(-1)
                    ,   a(0)
                    ,   b(0)
                    ,   result(false)
            {}
        };

#if 0
        template <class Container>
        void _check(const Container &new_pairs);
#endif

        void send_check_message(const pair<size_t,size_t> &pair);

        static void am_try_check_function(const AmContext * const context);
        static void am_complete_check_function(const AmContext * const context);
        static void am_complete_check_local_function(const AmContext * const context);

        void do_try_check_function(const TryCheck &try_check);
        void do_complete_check_function(const TryCheck &try_check);
        void do_complete_check_local_function(const TryCheck &try_check);

        SetPair pairs;
        int thread_rank;
        int check_id;
        AmHandle try_check;
        AmHandle complete_check;
        PthreadMutex server_mutex;
        Dispatcher<PthreadMutex> dispatcher;
        volatile int check_response;
        SetPair s_pairs_response;
        VecPair v_pairs_response;
};


#if 0
template <class Container>
void PairCheckGlobal::_check(const Container &new_pairs)
{
    typename Container::const_iterator;

    assert(check_response == 0);

    for (Container::const_iterator it=new_pairs.begin();
            it!=new_pairs.end(); ++it) {
        send_check_message(*it);
    }

    while (check_response > 0) {
#if !defined(THREADED)
        AmListenObjCodelet<NullMutex>* listenCodelet;
        if ((listenCodelet = theAm().amListeners[0]->progress()) != NULL) {
            listenCodelet->execute();
        }
#endif
    }
}
#endif

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKGLOBAL_H_ */
