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

void cntSort4Bkt(bucket_t *bkt, int bktSize);
void countSort(stnode_t *stNodes, int *srtNdIndex, int nStNodes, int maxSeqLen);

int bktCmp(const void *vp1, const void *vp2);
void printIntArr(int *srtNdIndex, int size);

#endif /* end of search.h */
