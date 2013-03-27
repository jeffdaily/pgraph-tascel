#ifndef BUCKET_H_
#define BUCKET_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "type.h"
#include "elib.h"


void initBkt(BKT *bkt, int bktSize);
int entryIndex(char *kmer, int k);
void slideWindow(char *str, int strLen, int sid, BKT *bkt, int bktSize, SUFFIX *sf, int k);
void buildBkt(SEQ *seqs, int nseqs, BKT *bkt, int bktSize, SUFFIX *sf, int sfSize, int k);


int printBkt(BKT *bkt, int bIndex);
int printBktList(SUFFIX *bktList);

#endif /* end of bucket.h */
