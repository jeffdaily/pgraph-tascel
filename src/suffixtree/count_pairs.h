#ifndef COUNT_PAIRS_H_
#define COUNT_PAIRS_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "stree.h"

unsigned long long count_pairs(stnode_t *stNodes, int *srtIndex, int nStNodes, int nSeqs, int EM, int *dup);

#endif /* end of pairs.h */
