/**
 * @file SuffixBucketsTascel.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SUFFIXBUCKETSTASCEL_H_
#define _PGRAPH_SUFFIXBUCKETSTASCEL_H_

#include <mpi.h>

#include <cstddef>
#include <vector>

#include <tascel.h>

#include "Suffix.hpp"
#include "SuffixBuckets.hpp"

using ::std::size_t;
using ::std::vector;

namespace pgraph {

class Parameters;
class SequenceDatabase;

class SuffixBucketsTascel : public SuffixBuckets
{
    public:
        /** @copydoc SuffixBuckets::SuffixBuckets */
        SuffixBucketsTascel(SequenceDatabase *sequences,
                      const Parameters &param,
                      MPI_Comm comm);

        virtual ~SuffixBucketsTascel();

        /** @copydoc SuffixBuckets::get */
        virtual Bucket* get(size_t bid);

        /** @copydoc SuffixBuckets::rem */
        virtual void rem(Bucket*);

        /** @copydoc SuffixBuckets::owns(size_t) */
        virtual bool owns(size_t bid) const;

        /** @copydoc SuffixBuckets::owns() */
        virtual const vector<size_t>& owns() const { return owned_buckets; }

        /** @copydoc SuffixBuckets::size_local() */
        virtual size_t size_local() const { return buckets_size; }

        /** @copydoc SuffixBuckets::size_nonempty() */
        virtual size_t size_nonempty() const { return n_nonempty; }

    private:
        /** For indexing into locally-stored Suffix allocation. */
        struct BucketMeta {
            size_t offset;
            size_t size;
            size_t bid;
        };

        /** Functor for sorting Suffix instances based on Bucket owner. */
        struct SuffixOwnerCompareFunctor {
            size_t comm_size;

            SuffixOwnerCompareFunctor(size_t comm_size)
                :   comm_size(comm_size)
            { }

            bool operator()(const Suffix &i, const Suffix &j) {
                size_t owner_i = i.bid % comm_size;
                size_t owner_j = j.bid % comm_size;
                return owner_i < owner_j;
            }
        };

        static bool SuffixBucketCompare(const Suffix &i, const Suffix &j) {
                return i.bid < j.bid;
        }

        ::tascel::AllocId aid_suffixes; /**< ID for suffixes allocation */
        ::tascel::AllocId aid_meta;     /**< ID for bucket meta allocation */
        Suffix *suffixes;               /**< array of all local suffixes */
        size_t suffixes_size;           /**< size of local suffixes array */
        BucketMeta *buckets;            /**< array of all local buckets */
        size_t buckets_size;            /**< size of local buckets array */
        size_t n_nonempty;              /**< global count of non-empty buckets */
        size_t count_remote_buckets;
        size_t count_remote_suffixes;
        vector<size_t> owned_buckets;
};


}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXBUCKETSTASCEL_H_ */
