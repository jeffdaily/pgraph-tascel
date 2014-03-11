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

#include "PairCheck.hpp"

using ::std::distance;
using ::std::inserter;
using ::std::set_difference;
using ::std::unique;

namespace pgraph {

class PairCheckSemiLocal : public PairCheck
{
    public:
        PairCheckSemiLocal() {}
        virtual ~PairCheckSemiLocal() {}

        virtual SetPair check(const SetPair &new_pairs) {
            SetPair ret_pairs;
            (void)set_difference(
                    new_pairs.begin(), new_pairs.end(),
                    s_pairs.begin(), s_pairs.end(),
                    inserter(ret_pairs, ret_pairs.end()));
            s_pairs.insert(ret_pairs.begin(), ret_pairs.end());
            return ret_pairs;
        }

        virtual VecPair check(const VecPair &new_pairs) {
            /* sort and unique-ify the incoming pairs */
            VecPair ret_pairs(new_pairs);
            sort(ret_pairs.begin(), ret_pairs.end());
            VecPair::iterator it_u = unique(ret_pairs.begin(), ret_pairs.end());
            ret_pairs.resize(distance(ret_pairs.begin(), it_u));
            /* set difference with already sorted local pairs */
            VecPair::iterator it_d = set_difference(
                    ret_pairs.begin(), ret_pairs.end(),
                    v_pairs.begin(), v_pairs.end(),
                    ret_pairs.begin());
            ret_pairs.resize(distance(ret_pairs.begin(), it_d));
            /* append, sort, and unique-ify new pairs */
            v_pairs.insert(v_pairs.end(), ret_pairs.begin(), ret_pairs.end());
            sort(v_pairs.begin(), v_pairs.end());
            it_u = unique(v_pairs.begin(), v_pairs.end());
            v_pairs.resize(distance(v_pairs.begin(), it_u));
            return ret_pairs;
        }

    private:
        SetPair s_pairs;
        VecPair v_pairs;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKSEMILOCAL_H_ */
