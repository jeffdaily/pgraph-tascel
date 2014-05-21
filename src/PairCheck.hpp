/**
 * @file PairCheck.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECK_H_
#define _PGRAPH_PAIRCHECK_H_

#include <cstddef>
#include <set>
#include <utility>
#include <vector>

#include "Pair.hpp"
#include "Stats.hpp"

using ::std::size_t;

namespace pgraph {

class PairCheck
{
    public:
        PairCheck() {}
        virtual ~PairCheck() {}

        virtual SetPair check(const SetPair &pairs) = 0;
        virtual VecPair check(const VecPair &pairs) = 0;
        virtual size_t size() = 0;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECK_H_ */
