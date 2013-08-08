/**
 * @file bucket.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_BUCKET_H_
#define _PGRAPH_BUCKET_H_

#include <cstddef>

#include "param.hpp"
#include "csequence.h"

namespace pgraph {

class SequenceDatabase;

/**
 * A suffix of a sequence_t, as a linked list item.
 */
typedef struct suffix {
    size_t sid;             /**< string id */
    size_t pid;             /**< position id */
    struct suffix *next;    /**< link to next suffix_t  */
} suffix_t;


/**
 * Bucket structure (linked list of related suffixes).
 */
typedef struct {
    suffix_t *suffixes; /**< linked list of zero or more suffixes */
    size_t size;        /**< number of suffixes */
} bucket_t;


/**
 * Suffix buckets for suffix tree processing.
 */
typedef struct {
    suffix_t *suffixes;     /**< array of all suffixes */
    size_t suffixes_size;   /**< size of suffixes array */
    bucket_t *buckets;      /**< array of all suffix buckets */
    size_t buckets_size;    /**< size of buckets array */
} suffix_buckets_t;


/**
 * Create and separate suffixes into appropriate suffix prefix buckets.
 *
 * @param[in] sequences to process
 * @param[in] param configuration parameters
 * @return suffixes and buckets
 */
suffix_buckets_t*
create_suffix_buckets_old(const sequences_t *sequences, param_t param);


/**
 * Free memory associated with suffixes and buckets.
 */
void free_suffix_buckets(suffix_buckets_t *suffix_buckets);


/**
 * Prints bucket buckets[index] to stdout.
 *
 * For debugging purposes.
 *
 * @param[in] buckets bucket pointers
 * @param[in] index for each bucket
 * @return the number of suffixes printed from the bucket
 */
size_t print_bucket(bucket_t *buckets, size_t index);


/**
 * Print a linked list of suffixes.
 *
 * @param[in] suffixes linked list of suffixes
 * @return the number of suffixes printed
 */
size_t print_suffixes(suffix_t *suffixes);

}; /* namespace pgraph */

#endif /* _PGRAPH_BUCKET_H_ */
