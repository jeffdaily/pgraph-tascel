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
#include <string>
#include <vector>

#include "Stats.hpp"

using ::std::size_t;
using ::std::string;
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
            ,   k(-1)
            ,   owner(-1)
        { }

        Suffix *suffixes;   /**< linked list of zero or more suffixes */
        size_t size;        /**< number of suffixes */
        size_t bid;         /**< bucket ID/index */
        int k;              /**< k-mer size */
        int owner;          /**< bucket owner */
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
         * @param[in] index of Bucket in remote Bucket array
         * @return the Bucket
         */
        virtual Bucket* get(int owner, size_t index) = 0;

        /**
         * Free the memory of the given Bucket instance.
         *
         * @param[in] bucket the Bucket
         */
        virtual void rem(Bucket *bucket) = 0;

        /**
         * Returns the number of buckets managed by this SuffixBuckets
         * instance.
         *
         * @return the size of this SuffixBuckets instance
         */
        virtual size_t size() const = 0;

        /**
         * Returns the number of buckets managed by this SuffixBuckets
         * instance on the current process.
         *
         * @return the number of process-local Bucket instances
         */
        virtual size_t size_local() const = 0;

        /**
         * Returns the number of non-NULL buckets managed by this
         * SuffixBuckets instance.
         *
         * @return the size of this SuffixBuckets instance
         */
        virtual size_t size_nonempty() const = 0;

        /**
         * Returns the (global) Stats for the bucket sizes.
         */
        virtual Stats stats_bucket_sizes() const = 0;

        /**
         * Returns Bucket index for the given Suffix string.
         *
         * Rather, the first 'k' characters of the given k-mer are the
         * prefix of the k-mer (suffix) where 'k' is the slide window
         * size as indicated by the stored Parameters instance. This
         * suffix prefix is calculated as a 0-based bucket index.
         *
         * @see bucket_kmer(size_t bid)
         * @param[in] kmer address of of k-mer (suffix) string
         * @return the Bucket index
         */
        size_t bucket_index(const char *kmer) const;

        /**
         * Returns Bucket index for the given Suffix string.
         *
         * Rather, the first 'k' characters of the given k-mer are the
         * prefix of the k-mer (suffix) where 'k' is the slide window
         * size as indicated by the stored Parameters instance. This
         * suffix prefix is calculated as a 0-based bucket index.
         *
         * @see bucket_kmer(size_t bid)
         * @param[in] kmer address of of k-mer (suffix) string
         * @param[in] k length of kmer
         * @return the Bucket index
         */
        size_t bucket_index(const char *kmer, int k) const;

        /**
         * Returns kmer for a given Bucket index.
         * @see bucket_index(const char*)
         * @param[in] bid bucket ID
         * @return string representing kmer
         */
        string bucket_kmer(size_t bid) const;

        /**
         * Returns kmer for a given Bucket index.
         * @see bucket_index(const char*)
         * @param[in] bid bucket ID
         * @param[in] k length of kmer
         * @return string representing kmer
         */
        string bucket_kmer(size_t bid, int k) const;
        
        static const size_t npos;

    protected:
        /**
         * Compare kmer against user-supplied prefix filters.
         *
         * @param[in] kmer address of of k-mer (suffix) string
         * @return true if the kmer should be skipped
         */
        bool filter_out(const char *kmer) const;

        /**
         * Compare kmer against user-supplied prefix filters.
         *
         * @param[in] kmer address of of k-mer (suffix) string
         * @param[in] k length of k-mer
         * @return true if the kmer should be skipped
         */
        bool filter_out(const char *kmer, int k) const;

        /**
         * power function for size_t data type
         */
        static size_t powz(size_t base, size_t n);

        MPI_Comm comm;                  /**< communicator */
        size_t comm_rank;               /**< communicator rank */
        size_t comm_size;               /**< communicator size */
        SequenceDatabase *sequences;    /**< sequences to process */
        const Parameters &param;        /**< configuration parameters */
        size_t SIGMA;                   /**< alphabet size */
        vector<size_t> alphabet_table;  /**< lookup table for alphabet index */
};

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXBUCKETS_H_ */
