/**
 * @file SuffixBuckets.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SUFFIXBUCKETS_H_
#define _PGRAPH_SUFFIXBUCKETS_H_

#include <mpi.h>

#include <cstddef>
#include <vector>

using ::std::size_t;
using ::std::vector;

namespace pgraph {

class Parameters;
class SequenceDatabase;
class Suffix;

/**
 * Bucket structure (linked list of related suffixes).
 */
class Bucket
{
    public:
        Bucket()
            :   suffixes(NULL)
            ,   size(0)
            ,   bid(0)
        { }

        Suffix *suffixes;   /**< linked list of zero or more suffixes */
        size_t size;        /**< number of suffixes */
        size_t bid;         /**< bucket ID/index */
};


class SuffixBuckets
{
    public:
        /**
         * Create and separate suffixes into appropriate suffix prefix buckets.
         * Distribut buckets evenly among available MPI ranks.
         *
         * @param[in] sequences to process
         * @param[in] param configuration parameters
         * @param[in] comm MPI communicator
         */
        SuffixBuckets(SequenceDatabase *sequences,
                      const Parameters &param,
                      MPI_Comm comm);

        /**
         * Free memory associated with suffixes and buckets.
         */
        virtual ~SuffixBuckets();

        /**
         * Retrieve a bucket, either remote or local.
         *
         * @param[in] bid Bucket index
         * @return the Bucket
         */
        virtual Bucket* get(size_t bid) = 0;

        /**
         * Free the memory of the given Bucket instance.
         *
         * @param[in] bucket the Bucket
         */
        virtual void rem(Bucket *bucket) = 0;

        /**
         * Returns true if the Bucket is owned by this process.
         *
         * That is, the Suffix instances that belong in the Bucket are
         * allocated on this process.
         *
         * @param[in] bid bucket index
         */
        virtual bool owns(size_t bid) const = 0;

        /**
         * Returns vector of all Bucket indexes owned by this process.
         *
         * That is, the Suffix instances that belong in the Bucket are
         * allocated on this process.
         *
         * @param[in] bid bucket index
         */
        virtual const vector<size_t> & owns() const = 0;

        /**
         * Returns the number of buckets managed by this SuffixBuckets
         * instance.
         *
         * @return the size of this SuffixBuckets instance
         */
        size_t size() const { return n_buckets; }

    protected:
        /**
         * Returns Bucket index for the given Suffix string.
         *
         * Rather, the first 'k' characters of the given k-mer are the
         * prefix of the k-mer (suffix) where 'k' is the slide window
         * size as indicated by the stored Parameters instance. This
         * suffix prefix is calculated as a 0-based bucket index.
         *
         * @param[in] kmer address of of k-mer (suffix) string
         * @return the Bucket index
         */
        size_t bucket_index(const char *kmer);

        /**
         * Compare kmer against user-supplied prefix filters.
         *
         * @param[in] kmer address of of k-mer (suffix) string
         * @return true if the kmer should be skipped
         */
        bool filter_out(const char *kmer);

        MPI_Comm comm;                  /**< communicator */
        size_t comm_rank;               /**< communicator rank */
        size_t comm_size;               /**< communicator size */
        size_t n_buckets;               /**< number of buckets */
        SequenceDatabase *sequences;    /**< sequences to process */
        const Parameters &param;        /**< configuration parameters */
        size_t SIGMA;                   /**< alphabet size */
        vector<size_t> alphabet_table;  /**< lookup table for alphabet index */
        
        static const size_t npos;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXBUCKETS_H_ */
