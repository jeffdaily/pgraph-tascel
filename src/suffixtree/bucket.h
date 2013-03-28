/**
 * @file bucket.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef BUCKET_H_
#define BUCKET_H_

#include <stddef.h>

#include "loadseq.h"


/**
 * Lookup table, used to represent suffix of every sequence.
 */
typedef struct suffix {
    int sid;                /**< string id */
    int pid;                /**< position id */
    struct suffix *next;    /**< ptr to next suffix  */
} suffix_t;


/**
 * Bucket structure.
 */
typedef struct {
    suffix_t *suffixes; /**< suffix array of zero or more suffixes */
    int size;           /**< number of suffixes */
} bucket_t;


/**
 * Bkt here is essentially counters. Init them to zero.
 *
 * @param[in] bucket array of buckets which contain suffixes
 * @param[in] size number of buckets
 */
void init_buckets(bucket_t *bucket, size_t size);


/**
 * Bucketing for sequence s <sid>.
 *
 * @param[in] s seqence itself
 * @param[in] s_len strlen(s) including ending '$'
 * @param[in] sid string id e.g. 0, 1, 2, ..., N-1
 * @param[out] buckets headers for each bucket; updated with new suffixes
 * @param[in] buckets_size bucket size
 * @param[out] suffixes suffix_t memory allocated for all buckets;
 *             updated with new suffixes used
 * @param[in] window_size slide window size
 */
void slide_window(const char *s, size_t s_len, int sid,
        bucket_t *buckets, size_t buckets_size, suffix_t *suffixes,
        int window_size);


/**
 * Buckets all sequences into buckets.
 *
 * @param[in] seqs fasta seqs
 * @param[in] nseqs #(fasta seqs)
 * @param[out] buckets global bucket shared by all
 * @param[in] buckets_size bucket size, which is <SIGMA^k>
 * @param[out] suffixes suffix_t memory allocated for all buckets to share
 * @param[in] suffixes_count how many suffixes have been used until now
 * @param[in] window_size slide window size
 */
void build_buckets(sequence_t *seqs, size_t nseqs,
        bucket_t *buckets, size_t buckets_size,
        suffix_t *suffixes, size_t suffix_count,
        int window_size);


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

#endif /* end of bucket.h */
