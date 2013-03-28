#include <assert.h>
#include <stdlib.h>

#include "bucket.h"
#include "dynamic.h" /* for SIGMA */
#include "loadseq.h"

/** @todo TODO doc */
static size_t _suffix_count = 0;

/**
 * This function can be optimized if higher performance is required.
 *
 * @param[in] kmer address of of k-mer string
 * @param[in] k slide window size
 * @return TODO todo
 */
static size_t entry_index(const char *kmer, int k)
{
    int i;
    size_t value = 0;

    for (i=0; i<k; ++i) {
        value = value*SIGMA + (kmer[i] - 'A');
    }

    return value;
}


void init_buckets(bucket_t *buckets, size_t size)
{
    size_t i;
    for(i=0; i<size; ++i){
        buckets[i].suffixes = NULL; 
        buckets[i].size = 0; 
    }
}


void slide_window(const char *s, size_t s_len, int sid,
        bucket_t *buckets, size_t buckets_size, suffix_t *suffixes,
        int window_size)
{
    int i;
    int stop_index = s_len - window_size - 1;
   
    for(i=0; i<=stop_index; ++i){
        size_t bucket_index = entry_index(s+i, window_size);
        assert(bucket_index < buckets_size);
        /* prefixed in the bucket list */
        suffixes[_suffix_count].sid = sid;
        suffixes[_suffix_count].pid = i;
        suffixes[_suffix_count].next = buckets[bucket_index].suffixes;
        buckets[bucket_index].suffixes = &suffixes[_suffix_count];
        buckets[bucket_index].size++;
        _suffix_count++;
    }
}


void build_buckets(sequence_t *seqs, size_t nseqs,
        bucket_t *buckets, size_t buckets_size,
        suffix_t *suffixes, size_t suffixes_count,
        int window_size)
{
    size_t i;
    init_buckets(buckets, buckets_size);
    
    /* slide k-mers for every seqs, and bucket them */
    for(i=0; i<nseqs; i++){
        slide_window(seqs[i].str, seqs[i].strLen, i,
                buckets, buckets_size, suffixes, window_size);
    }

    printf("_suffix_count=%zu\n", _suffix_count);
    assert(_suffix_count == suffixes_count);
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
    for(p = suffixes; p != NULL; p = p->next){
        printf("[%d, %d]\t", p->sid, p->pid);
        ++i;
    }
    printf("\n");

    return i;
}
