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
        virtual Bucket* get(int owner, size_t index);

        /** @copydoc SuffixBuckets::rem */
        virtual void rem(Bucket*);

        /** @copydoc SuffixBuckets::size() */
        virtual size_t size() const { return buckets_size_global; }

        /** @copydoc SuffixBuckets::size_local() */
        virtual size_t size_local() const { return buckets_size; }

        /** @copydoc SuffixBuckets::size_nonempty() */
        virtual size_t size_nonempty() const { return n_nonempty; }

        /** @copydoc SuffixBuckets::stats_bucket_sizes() */
        virtual Stats stats_bucket_sizes() const { return buckets_stats; }

    private:
        /** For indexing into locally-stored Suffix allocation. */
        struct BucketMeta {
            size_t offset;
            size_t size;
            size_t bid;
            int k;

            BucketMeta() : offset(0U), size(0U), bid(0U), k(-1) {}
            BucketMeta(size_t offset, size_t size, size_t bid, int k)
                : offset(offset), size(size), bid(bid), k(k) {}
        };

        void refine_bucket(BucketMeta &bucket, vector<BucketMeta> &buckets);

        ::tascel::AllocId aid_suffixes; /**< ID for suffixes allocation */
        ::tascel::AllocId aid_meta;     /**< ID for bucket meta allocation */
        Suffix *suffixes;               /**< array of all local suffixes */
        size_t suffixes_size;           /**< size of local suffixes array */
        BucketMeta *buckets;            /**< array of all local buckets */
        size_t buckets_size;            /**< size of local buckets array */
        size_t buckets_size_global;     /**< size of global buckets array */
        size_t n_nonempty;              /**< global count of non-empty buckets */
        Stats buckets_stats;            /**< stats for bucket sizes */
        size_t count_remote_buckets;
        size_t count_remote_suffixes;
};


}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXBUCKETSTASCEL_H_ */
