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
        virtual VecPair check(const VecPair &pairs);

    //protected:
    public:
        bool send_check_message(const pair<size_t,size_t> &pair);
        void bulk_send_check_message(const SetPair &pairs, SetPair &result);

        static void try_check_function(const AmContext * const context);
        static void check_complete_function(const AmContext * const context);
        static void check_complete_local_function(const AmContext * const context);

        static void bulk_try_check_function(const AmContext * const context);
        static void bulk_check_complete_function(const AmContext * const context);
        static void bulk_check_complete_local_function(const AmContext * const context);

        PthreadMutex mutex;
        int thd;
        AllocId alloc_id;
        AmHandle try_check;
        AmHandle check_complete;
        AmHandle bulk_try_check;
        AmHandle bulk_check_complete;
        SetPair s_pairs;
        VecPair v_pairs;
        Dispatcher<NullMutex> dispatcher;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKGLOBAL_H_ */
