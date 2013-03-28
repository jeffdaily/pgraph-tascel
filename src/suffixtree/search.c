/**
 * @file search.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "bucket.h"
#include "elib.h"
#include "search.h"
#include "stree.h"


void count_sort_buckets(bucket_t *buckets, size_t n_buckets)
{
    int *count = NULL;
    size_t i = 0;
    bucket_t *new_buckets = NULL;
    size_t max_bucket_size = 0;

    /* find the maximum size for n_buckets */
    for (i=0; i<n_buckets; ++i) {
        if (buckets[i].size > max_bucket_size) {
            max_bucket_size = buckets[i].size;
        }
    }
    /* +1 for counters *count */
    ++max_bucket_size;

    count = emalloc(max_bucket_size*sizeof(int));
    new_buckets = emalloc(n_buckets*sizeof(bucket_t));

    /* init counter */
    for (i=0; i<max_bucket_size; ++i){
        count[i] = 0;
    }
    
    /* counting */
    for (i=0; i<n_buckets; ++i) {
        if (buckets[i].size >= max_bucket_size) {
            printf("[Cnt=%zu, mSize=%zu]\n", buckets[i].size, max_bucket_size); 
            exit(0);
        }
        ++count[buckets[i].size];
    }

    /* prefix sum in a reverse way */
    for (i=max_bucket_size-2; i>0; --i) {
        count[i] += count[i+1];
    }
    count[i] += count[i+1]; /* for the i == 0 case */
    
    /* put the buckets[] into new_buckets[] in the descending order */
    for (i=0; i<n_buckets; ++i) {
        new_buckets[count[buckets[i].size]-1] = buckets[i];
        --count[buckets[i].size];
    }

    /* copy descending order back */
    for (i=0; i<n_buckets; ++i){
        buckets[i] = new_buckets[i];
    }
    
    free(count); 
    free(new_buckets); 
}


/* -----------------------------------------------------------*
 * counting sort based on the depth of the stNodes[i].
 * 
 * @param stNodes -
 * @param srtIndex -
 * @param nStNodes -
 * @param maxSeqLen -
 * -----------------------------------------------------------*/
void countSort(stnode_t *stNodes, int *srtIndex, int nStNodes, int maxSeqLen){
    int i;
    int *cnt = NULL;
    int depth;

    cnt = emalloc(maxSeqLen*(sizeof *cnt));

    /* init counters */
    for(i = 0; i < maxSeqLen; i++){
        cnt[i] = 0;
    } 

    /* fill in counters */
    for(i = 0; i < nStNodes; i++){
        depth = stNodes[i].depth; /* depth starts from 0 */
        assert(depth < maxSeqLen);
        cnt[depth]++;
    }
    
    /* suffix sum to sort in a reverse way */
    for(i = maxSeqLen-2; i >= 0; i--){
       cnt[i] += cnt[i+1]; 
    }
    
    /* store the sorted index into srtIndex */
    for(i = 0; i < nStNodes; i++){
        depth = stNodes[i].depth;
        srtIndex[cnt[depth]-1] = i;
        cnt[depth]--;
    }

    free(cnt);
}


/* -----------------------------------------------------------*
 * comparison function for qsort().
 * NOTE: it is sorted in a descending way.
 * 
 * @param vp1 -
 * @param vp2 -
 * -----------------------------------------------------------*/
int bktCmp(const void *vp1, const void *vp2){
    bucket_t *b1 = (bucket_t *)vp1;
    bucket_t *b2 = (bucket_t *)vp2;

    if(b1->size > b2->size){
        return -1;
    }else if (b1->size == b2->size){
        return 0;
    }else{
        return 1;
    }
}


void printIntArr(int *srtIndex, int size){
    int i;
    for(i = 0; i < size; i++){
        printf("strNdIndex[%d]=%d\n", i, srtIndex[i]);
    }
}
