/**
 * @file search.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEARCH_H_
#define SEARCH_H_

#include "bucket.h"
#include "stree.h"


/**
 * Sorts the buckets according to their size. 
 *
 * @param[in] buckets array of all buckets to be sorted
 * @param[in] n_buckets #(buckets)
 */
void count_sort_buckets(bucket_t *buckets, size_t n_buckets);


void countSort(stnode_t *stNodes, int *srtNdIndex, int nStNodes, int maxSeqLen);

int bktCmp(const void *vp1, const void *vp2);
void printIntArr(int *srtNdIndex, int size);

#endif /* end of search.h */
