/**
 * @file PairCheckLocal.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2014 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PAIRCHECKLOCAL_H_
#define _PGRAPH_PAIRCHECKLOCAL_H_

#include "PairCheck.hpp"

namespace pgraph {

class PairCheckLocal : public PairCheck
{
    public:
        PairCheckLocal() {}
        virtual ~PairCheckLocal() {}

        virtual SetPair check(const SetPair &pairs) {
            return pairs;
        }
        virtual VecPair check(const VecPair &pairs) {
            return pairs;
        }
        virtual size_t size() {
            return 0;
        }
};

}; /* namespace pgraph */

#endif /* _PGRAPH_PAIRCHECKLOCAL_H_ */
