/**
 * @file Suffix.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SUFFIX_H_
#define _PGRAPH_SUFFIX_H_

#include <ostream>

namespace pgraph {

/**
 * A suffix of a Sequence.
 */
class Suffix
{
    public:
        Suffix()
            :   sid(0)
            ,   pid(0)
            ,   bid(0)
            ,   next(NULL)
        { }

        Suffix(size_t suffix_id,
                size_t offset,
                size_t bucket_id,
                Suffix *next=NULL)
            :   sid(suffix_id)
            ,   pid(offset)
            ,   bid(bucket_id)
            ,   next(next)
        { }

        size_t sid;     /**< suffix id */
        size_t pid;     /**< position id i.e. offset */
        size_t bid;     /**< bucket id */
        Suffix *next;   /**< link to next Suffix  */
};

inline ::std::ostream& operator<<(::std::ostream &os, const Suffix &suffix)
{
    return os << "(" << suffix.sid
        << "," << suffix.pid
        << "," << suffix.bid
        << ")";
}

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIX_H_ */
