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
                      bool dumb_split=false);

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
        size_t first_bucket;
        size_t last_bucket;
        std::vector<size_t> bucket_size;
        std::vector<size_t> bucket_owner;
        std::vector<size_t> bucket_offset;
        std::vector<Suffix*> bucket_address;
        size_t bucket_size_total;
        PthreadMutex mutex;
        bool dumb_split;
};

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXBUCKETS_H_ */
