/**
 * @file PairCheckSemiLocal.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECKSEMILOCAL_H_
#define _PGRAPH_PAIRCHECKSEMILOCAL_H_

#include <algorithm>
#include <iterator>

#include <mpi.h> /* for MPI_Wtime() */

#include "PairCheck.hpp"

using ::std::inserter;
using ::std::set_difference;

namespace pgraph {

class PairCheckSemiLocal : public PairCheck
{
    public:
        PairCheckSemiLocal() {}
        virtual ~PairCheckSemiLocal() {}

        virtual SetPair check(const SetPair &new_pairs) {
            double t = MPI_Wtime();
            SetPair ret_pairs;
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
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKSEMILOCAL_H_ */
