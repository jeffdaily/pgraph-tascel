/**
 * @file SuffixBuckets.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SUFFIXBUCKETS_H_
#define _PGRAPH_SUFFIXBUCKETS_H_

#include <mpi.h>

#include <cstddef>
#include <iostream>
#include <vector>

#include <tascel.h>

#include "Parameters.hpp"

using tascel::PthreadMutex;

namespace pgraph {

class Parameters;
class SequenceDatabase;

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

        Suffix(size_t suffix_id, size_t offset, size_t bucket_id, Suffix *next=NULL)
            :   sid(suffix_id)
            ,   pid(offset)
            ,   bid(0)
            ,   next(next)
        { }

        size_t sid;     /**< suffix id */
        size_t pid;     /**< position id i.e. offset */
        size_t bid;     /**< bucket id */
        Suffix *next;   /**< link to next Suffix  */
};

ostream& operator<<(ostream &os, const Suffix &suffix);



/**
 * Bucket structure (linked list of related suffixes).
 */
class Bucket
{
    public:
        Bucket()
            :   suffixes(NULL)
            ,   size(0)
        { }

        Suffix *suffixes;   /**< linked list of zero or more suffixes */
        size_t size;        /**< number of suffixes */
};

enum SplitEnum {SPLIT_DUMB=0, SPLIT_SUFFIXES=1, SPLIT_SORTED=2};

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
                      MPI_Comm comm,
                      const SplitEnum &split_type=SPLIT_SORTED);

        /**
         * Free memory associated with suffixes and buckets.
         */
        ~SuffixBuckets();

        Suffix *get(size_t bid);

        bool owns(size_t bid) const;

    //private:
        int comm_rank;
        int comm_size;
        SequenceDatabase *sequences;    /**< sequences to process */
        Parameters param;       /**< configuration parameters */
        Suffix *suffixes;       /**< array of all suffixes */
        size_t suffixes_size;   /**< size of suffixes array */
        Bucket *buckets;        /**< array of all suffix buckets */
        size_t buckets_size;    /**< size of buckets array */
        std::vector<size_t> my_buckets;
        std::vector<size_t> bucket_size;
        std::vector<size_t> bucket_owner;
        std::vector<size_t> bucket_offset;
        Suffix **bucket_address;
        size_t bucket_size_total;
        PthreadMutex mutex;
        SplitEnum split_type;

        bool operator()(const Suffix &i, const Suffix &j) {
            return bucket_owner[i.bid] < bucket_owner[j.bid];
        }
};


class Bucket2
{
    public:
        Bucket2()
            :   offset(0)
            ,   size(0)
        { }

        size_t offset;
        size_t size;
};


class SuffixBuckets2
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
        SuffixBuckets2(SequenceDatabase *sequences,
                      const Parameters &param,
                      MPI_Comm comm);

        /**
         * Free memory associated with suffixes and buckets.
         */
        ~SuffixBuckets2();

        Bucket* get(size_t bid);

        void rem(size_t bid, Bucket*);

        bool owns(size_t bid) const;

    //private:
        int comm_rank;
        int comm_size;
        size_t n_buckets;
        SequenceDatabase *sequences;    /**< sequences to process */
        Parameters param;               /**< configuration parameters */
        Suffix **suffixes_remote;       /**< addresses of remote suffixes */
        Suffix *suffixes;               /**< array of all local suffixes */
        size_t suffixes_size;           /**< size of local suffixes array */
        Bucket2 **buckets_remote;       /**< addresses of remote buckets */
        Bucket2 *buckets;               /**< array of all local buckets */
        size_t buckets_size;            /**< size of local buckets array */
        PthreadMutex mutex;
};


}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXBUCKETS_H_ */
