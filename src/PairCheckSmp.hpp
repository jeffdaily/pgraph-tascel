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

#include <mpi.h> /* for MPI_Wtime() */

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
            double t = MPI_Wtime();
            SetPair ret_pairs;
            LockGuard<PthreadMutex> guard(mutex);
            (void)set_difference(
                    new_pairs.begin(), new_pairs.end(),
                    pairs.begin(), pairs.end(),
                    inserter(ret_pairs, ret_pairs.end()));
            pairs.insert(ret_pairs.begin(), ret_pairs.end());
            time.push_back(MPI_Wtime() - t);
            return ret_pairs;
        }

        virtual VecPair check(const VecPair &new_pairs) {
            double t = MPI_Wtime();
            VecPair ret_pairs(new_pairs.size());
            LockGuard<PthreadMutex> guard(mutex);
            VecPair::iterator it = set_difference(
                    new_pairs.begin(), new_pairs.end(),
                    pairs.begin(), pairs.end(),
                    ret_pairs.begin());
            ret_pairs.resize(it-ret_pairs.begin());
            pairs.insert(ret_pairs.begin(), ret_pairs.end());
            time.push_back(MPI_Wtime() - t);
            return ret_pairs;
        }

    private:
        SetPair pairs;
        PthreadMutex mutex;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKSMP_H_ */
