/**
 * @file bucket.cpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "bucket.hpp"
#include "constants.h"
#include "csequence.h"

namespace pgraph {

/** power function for size_t */
size_t powz(size_t base, size_t n)
{
    size_t p = 1;

    for(/*empty*/; n > 0; --n) {
        assert(p < SIZE_MAX/base);
        p *= base;
    }

    return p;
}


/**
 * This function can be optimized if higher performance is required.
 *
 * @param[in] kmer address of of k-mer string
 * @param[in] k slide window size
 * @return TODO todo
 */
static inline size_t entry_index(const char *kmer, int k)
{
    int i;
    size_t value = 0;

    for (i = 0; i < k; ++i) {
        value = value * SIGMA + (kmer[i] - 'A');
    }

    return value;
}


suffix_buckets_t*
create_suffix_buckets_old(const sequences_t *sequences, param_t param)
{
    suffix_buckets_t *suffix_buckets = NULL;
    suffix_t *suffixes = NULL;
    bucket_t *buckets = NULL;
    size_t n_suffixes = 0;
    size_t n_buckets = 0;
    size_t i = 0;
    size_t suffix_index = 0;

    /* allocate buckets */
    n_buckets = powz(SIGMA, param.window_size);
    buckets = new bucket_t[n_buckets];
    if (NULL == buckets) {
        perror("create_suffix_buckets: malloc buckets");
        exit(EXIT_FAILURE);
    }

    /* initialize buckets */
    for (i = 0; i < n_buckets; ++i) {
        buckets[i].suffixes = NULL;
        buckets[i].size = 0;
    }

    /* allocate suffixes */
    n_suffixes = sequences->n_chars - sequences->size * param.window_size;
    suffixes = new suffix_t[n_suffixes];
    if (NULL == suffixes) {
        perror("create_suffix_buckets: malloc suffixes");
        exit(EXIT_FAILURE);
    }

    /* initialize suffixes */
    for (i = 0; i < n_suffixes; ++i) {
        suffixes[i].sid = 0;
        suffixes[i].pid = 0;
        suffixes[i].next = NULL;
    }

    /* slide k-mers for every sequence and bucket them */
    for (i = 0; i < sequences->size; ++i) {
        sequence_t *sequence = NULL;
        size_t stop_index = 0;
        size_t j = 0;

        sequence = &(sequences->seq[i]);
        stop_index = sequence->size - param.window_size - 1;

        assert(sequence->size >= param.window_size);

        for (j = 0; j <= stop_index; ++j) {
            size_t bucket_index = entry_index(
                    sequence->str + j, param.window_size);
            assert(bucket_index < n_buckets);
            /* prefixed in the suffix list for the given bucket */
            suffixes[suffix_index].sid = i;
            suffixes[suffix_index].pid = j;
            suffixes[suffix_index].next = buckets[bucket_index].suffixes;
            buckets[bucket_index].suffixes = &suffixes[suffix_index];
            buckets[bucket_index].size++;
            suffix_index++;
        }
    }

    /* return values */
    suffix_buckets = new suffix_buckets_t;
    if (NULL == suffix_buckets) {
        perror("create_suffix_buckets: malloc suffix_buckets");
        exit(EXIT_FAILURE);
    }
    suffix_buckets->suffixes = suffixes;
    suffix_buckets->suffixes_size = n_suffixes;
    suffix_buckets->buckets = buckets;
    suffix_buckets->buckets_size = n_buckets;

    return suffix_buckets;
}


void free_suffix_buckets(suffix_buckets_t *suffix_buckets)
{
    delete [] suffix_buckets->suffixes;
    delete [] suffix_buckets->buckets;
    delete suffix_buckets;
}


size_t print_bucket(bucket_t *buckets, size_t index)
{
    return print_suffixes(buckets[index].suffixes);
}


size_t print_suffixes(suffix_t *suffixes)
{
    suffix_t *p = NULL;
    size_t i = 0;

    printf("->");
    for (p = suffixes; p != NULL; p = p->next) {
        printf("[%zu, %zu]\t", p->sid, p->pid);
        ++i;
    }
    printf("\n");

    return i;
}

}; /* namespace pgraph */

