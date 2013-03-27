#ifndef STREE_H_
#define STREE_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "type.h"
#include "lib.h"
#include "elib.h"
#include "bucket.h"
#include "loadseq.h"
#include "search.h"
#include "pairs.h"
#include "uFind.h"

enum {ERROR = -1, DOL_END = -2};


int buildTree(SEQ *seqs, int nseqs, SUFFIX *bktList, int depth, int k, STNODE *stNodes, int *stIndex);
void compLset(SUFFIX *bktList, SEQ *seqs, int nseqs, SUFFIX **lset);
void procBkt(SEQ *seqs, int nseqs, SUFFIX *bktList, int blSize, int maxSeqLen, UF *uSet, PARAM *param, int ind);
void buildForest(BKT *bkt, int bktSize, SEQ *seqs, int nseqs, int maxSeqLen, UF *uSet, PARAM *param);
int nextDiffPos(SEQ *seqs, SUFFIX *lset, int depth, int k);


int lsetSize(SUFFIX *lset);
void freeLset(SUFFIX *lset);
void freeStNodes(STNODE *stNodes, int size);


void printStnodes(STNODE *stNodes, int stIndex, int blSize, int ind);
void printLset(SUFFIX **lset);
void printCnt(int *cnt, int size);
void printStNode(STNODE *stNodes, int i);
void printStNodes(STNODE *stNodes, int offset);

#endif /* end of stree.h */
