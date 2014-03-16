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
            for (SetPair::const_iterator it=new_pairs.begin();
                    it!=new_pairs.end(); ++it) {
                pair<SetPair::iterator,bool> result = s_pairs.insert(*it);
                if (result.second) {
                    (void)ret_pairs.insert(*it);
                }
            }
            return ret_pairs;
        }

        virtual VecPair check(const VecPair &new_pairs) {
            LockGuard<PthreadMutex> guard(mutex);
            VecPair ret_pairs;
            for (VecPair::const_iterator it=new_pairs.begin();
                    it!=new_pairs.end(); ++it) {
                pair<SetPair::iterator,bool> result = s_pairs.insert(*it);
                if (result.second) {
                    (void)ret_pairs.push_back(*it);
                }
            }
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
