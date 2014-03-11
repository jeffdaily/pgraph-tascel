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

using ::tascel::AllocId;
using ::tascel::AmArg;
using ::tascel::AmContext;
using ::tascel::AmHandle;
using ::tascel::RmaPtr;
using ::tascel::PthreadMutex;

namespace pgraph {

class PairCheckGlobal : public PairCheck
{
    public:
        PairCheckGlobal();
        virtual ~PairCheckGlobal();

        virtual SetPair check(const SetPair &pairs);
        virtual VecPair check(const VecPair &pairs);

    protected:
        struct PairCheckArg : public AmArg {
            const RmaPtr ptr;
            PairCheckArg(const RmaPtr& p) : AmArg(), ptr(p) {}
        };

        bool send_check_message(const pair<size_t,size_t> &pair);
        bool do_check(size_t s_pair[2]);

        static void amLocalClient(const AmContext * const context);
        static void amPostPutServer(const AmContext * const context);
        static void amRemoteClient(const AmContext * const context);
        static void amLocalServer(const AmContext * const context);
        static void amPrePutServer(const AmContext * const context);

        PthreadMutex mutex;
        AllocId alloc_id;
        AmHandle am_handle;
        SetPair s_pairs;
        VecPair v_pairs;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKGLOBAL_H_ */
