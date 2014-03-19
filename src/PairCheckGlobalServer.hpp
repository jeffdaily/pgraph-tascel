/**
 * @file PairCheckGlobalServer.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECKGLOBALSERVER_H_
#define _PGRAPH_PAIRCHECKGLOBALSERVER_H_

#include <tascel.h>

#include "PairCheck.hpp"

using ::tascel::AmContext;
using ::tascel::AmHandle;
using ::tascel::Dispatcher;
using ::tascel::PthreadMutex;
using ::tascel::NullMutex;

namespace pgraph {

class PairCheckGlobalServer : public PairCheck
{
    public:
        PairCheckGlobalServer(int thd);
        virtual ~PairCheckGlobalServer();

        virtual SetPair check(const SetPair &pairs);
        virtual VecPair check(const VecPair &pairs);

        virtual size_t size();

    //protected:
    public:
        void bulk_send_check_message(const SetPair &pairs, SetPair &result);

        static void bulk_try_check_function(const AmContext * const context);
        static void bulk_try_check_local_function(const AmContext * const context);
        static void bulk_check_complete_function(const AmContext * const context);
        static void bulk_check_complete_local_function(const AmContext * const context);

        PthreadMutex mutex;
        int thd;
        AmHandle bulk_try_check;
        AmHandle bulk_check_complete;
        Dispatcher<NullMutex> dispatcher;
        size_t last_size;
        volatile int server_response;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKGLOBALSERVER_H_ */
