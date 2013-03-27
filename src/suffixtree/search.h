#ifndef SEARCH_H_
#define SEARCH_H_

#include <stdio.h>
#include <assert.h>
#include "type.h"
#include "elib.h"

void cntSort4Bkt(BKT *bkt, int bktSize);
void countSort(STNODE *stNodes, int *srtNdIndex, int nStNodes, int maxSeqLen);

int bktCmp(const void *vp1, const void *vp2);
void printIntArr(int *srtNdIndex, int size);

#endif /* end of search.h */
