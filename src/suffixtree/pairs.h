#ifndef PAIRS_H_
#define PAIRS_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "type.h"
#include "stree.h"
#include "dynamic.h"
#include "uFind.h"

void printPairs();
void genPairs(STNODE *stNodes, int *srtIndex, int nStNodes, SEQ *seqs, int nSeqs,
                int maxSeqLen, UF *uSet, int *dup, PARAM *param);

void procLeaf(SUFFIX **lset, SEQ *seqs, int nSeqs, CELL **tbl, int **ins, int **del, UF *uSet, PARAM *param);

int isEdge(CELL *result, char *s1, int s1Len, char *s2, int s2Len, PARAM *param);

#endif /* end of pairs.h */
