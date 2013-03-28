#ifndef PAIRS_H_
#define PAIRS_H_

#include "dynamic.h"
#include "loadseq.h"
#include "stree.h"

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
                int maxSeqLen, ufind_t *uSet, int *dup, param_t *param);

void procLeaf(suffix_t **lset, sequence_t *seqs, int nSeqs, cell_t **tbl, int **ins, int **del, ufind_t *uSet, param_t *param);

int isEdge(cell_t *result, char *s1, int s1Len, char *s2, int s2Len, param_t *param);

#endif /* end of pairs.h */
