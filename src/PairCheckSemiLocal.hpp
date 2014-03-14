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

        virtual size_t size() {
            return s_pairs.size();
        }

    private:
        SetPair s_pairs;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKSEMILOCAL_H_ */
