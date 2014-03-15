/**
 * @file PairCheckGlobal.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECKGLOBAL_H_
#define _PGRAPH_PAIRCHECKGLOBAL_H_

#include <vector>

#include <tascel.h>

#include "PairCheck.hpp"

using ::std::vector;
using ::tascel::AllocId;
using ::tascel::AmArg;
using ::tascel::AmContext;
using ::tascel::AmHandle;
using ::tascel::Dispatcher;
using ::tascel::PthreadMutex;
using ::tascel::NullMutex;

namespace pgraph {

class PairCheckGlobal : public PairCheck
{
    public:
        PairCheckGlobal(int thd);
        virtual ~PairCheckGlobal();

        virtual SetPair check(const SetPair &pairs);

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
        SetPair s_pairs;
        Dispatcher<NullMutex> dispatcher;
        volatile int server_response;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKGLOBAL_H_ */
