/**
 * @file PairCheckSmp.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECKSMP_H_
#define _PGRAPH_PAIRCHECKSMP_H_

#include <algorithm>
#include <iterator>

#include <tascel.h>

#include "PairCheck.hpp"

using ::std::inserter;
using ::std::set_difference;
using ::tascel::LockGuard;
using ::tascel::PthreadMutex;

namespace pgraph {

class PairCheckSmp : public PairCheck
{
    public:
        PairCheckSmp() {}
        virtual ~PairCheckSmp() {}

        virtual SetPair check(const SetPair &new_pairs) {
            LockGuard<PthreadMutex> guard(mutex);
            SetPair ret_pairs;
            (void)set_difference(
                    new_pairs.begin(), new_pairs.end(),
                    s_pairs.begin(), s_pairs.end(),
                    inserter(ret_pairs, ret_pairs.end()));
            s_pairs.insert(ret_pairs.begin(), ret_pairs.end());
            return ret_pairs;
        }

        virtual size_t size() {
            LockGuard<PthreadMutex> guard(mutex);
            return s_pairs.size();
        }

    private:
        SetPair s_pairs;
        PthreadMutex mutex;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKSMP_H_ */
