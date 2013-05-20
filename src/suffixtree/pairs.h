/**
 * @file pairs.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef PAIRS_H_
#define PAIRS_H_

#include "dynamic.h"
#include "loadseq.h"
#include "stree.h"
//#include "uFind.h"

/**
 * Prints pair statistics to stdout.
 */
void print_pairs();


/**
 * Generate promising pairs for alignment.
 *
 * @param[in] stNodes -
 * @param[in] srtIndex -
 * @param[in] nStNodes -
 * @param[in] seqs -
 * @param[in] nSeqs -
 * @param[in] maxSeqLen -
 * @param[in] uSet -
 * @param[in] dup -
 * @param[in] param -
 */
void genPairs(stnode_t *stNodes, int *srtIndex, int nStNodes, sequence_t *seqs, int nSeqs,
              int maxSeqLen, int *dup, param_t *param);

/**
 * This function implements the pair generation algorithm for leaf nodes.
 *
 *    BEGIN - intra/inter cross. O/W - intra cross.
 *
 * @param lset -
 * @param seqs -
 * @param nSeqs -
 * @param tbl -
 * @param ins -
 * @param del -
 */
void procLeaf(suffix_t **lset, sequence_t *seqs, int nSeqs, cell_t **tbl, int **ins, int **del, param_t *param, int *dup);

int isEdge(cell_t *result, char *s1, int s1Len, char *s2, int s2Len, param_t *param);

#endif /* end of pairs.h */
