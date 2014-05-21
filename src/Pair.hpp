/**
 * @file Pair.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIR_H_
#define _PGRAPH_PAIR_H_

#include <cstddef>
#include <set>
#include <utility>
#include <vector>

using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::vector;

namespace pgraph {

typedef pair<size_t,size_t> Pair;
typedef set<Pair> SetPair;
typedef vector<Pair> VecPair;

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIR_H_ */
