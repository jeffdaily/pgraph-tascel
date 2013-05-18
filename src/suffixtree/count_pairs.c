/*
 * $Rev: 76 $
 * $Date: 2011-05-11 16:00:44 -0700 (Wed, 11 May 2011) $
 * $Author: andy.cj.wu@gmail.com $
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * ----------------------------------------------------------------
 *
 */

#include "count_pairs.h"

/*---------------------------------------------------------------*
 * This function implements the pair generation algorithm for leaf
 * nodes.
 *
 *    BEGIN - intra/inter cross. O/W - intra cross.
 *
 * @param lset -
 *---------------------------------------------------------------*/
static void procLeaf(suffix_t **lset, unsigned long long *nPairs)
{
    int i;
    int j;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int f1, f2;

    for (i = 0; i < SIGMA; i++) {
        if (lset[i]) {
            if (i == BEGIN - 'A') { /* inter cross */
                for (p = lset[i]; p != NULL; p = p->next) {
                    for (q = p->next; q != NULL; q = q->next) {
                        f1 = p->sid;
                        f2 = q->sid;

                        if (f1 == f2) {
                            continue;
                        }

                        (*nPairs)++;
                    }
                }
            }

            /* intra cross */
            for (j = i + 1; j < SIGMA; j++) {
                if (lset[j]) {
                    for (p = lset[i]; p != NULL; p = p->next) {
                        for (q = lset[j]; q != NULL; q = q->next) {
                            f1 = p->sid;
                            f2 = q->sid;

                            if (f1 == f2) {
                                continue;
                            }

                            (*nPairs)++;
                        }
                    }
                }
            }
        }
    }
}

/*---------------------------------------------------------------*
 * Generate promising pairs, and enBuf them into pBuf[]
 *
 * @param stNodes -
 * @param srtIndex -
 * @param nStNodes -
 * @param nSeqs -
 * @param EM - exact match length cutoff
 * @param dup - duplication arrary for redundant pairs checking
 *
 *---------------------------------------------------------------*/
unsigned long long count_pairs(stnode_t *stNodes, int *srtIndex, int nStNodes, int nSeqs,
                               int EM, int *dup)
{
    int i;
    int j;
    int r;
    stnode_t *stnode = NULL;
    int sIndex, eIndex;
    int m, n, s, t;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int f1, f2;
    unsigned long long nPairs = 0;

    /* srtIndex maintain an order of NON-increasing depth of stNodes[] */
    for (i = 0; i < nStNodes; i++) {
        sIndex = srtIndex[i];
        stnode = &stNodes[sIndex];

#ifdef DEBUG
        printf("stNode->depth=%d, stnode->rLeaf=%ld, sIndex=%ld\n", stnode->depth, stnode->rLeaf, sIndex);
#endif

        if (stnode->depth >= EM - 1) {
            if (stnode->rLeaf == sIndex) { /* leaf node */
                procLeaf(stnode->lset, &nPairs);
            }
            else {                       /* internal node */
                eIndex = stnode->rLeaf;

                /* init dup[] for the internal node */
                for (r = 0; r < nSeqs; r++) {
                    dup[r] = -1;
                }

                /* pairs generation loop for internal node */
                for (m = sIndex + 1; m < eIndex; m = stNodes[m].rLeaf + 1) {
                    for (n = stNodes[m].rLeaf + 1; n <= eIndex; n = stNodes[n].rLeaf + 1) {
                        for (s = 0; s < SIGMA; s++) {
                            if (stNodes[m].lset[s]) {
                                for (t = 0; t < SIGMA; t++) {
                                    if (stNodes[n].lset[t]) {
                                        if (s != t) {
                                            for (p = stNodes[m].lset[s]; p != NULL; p = p->next) {

                                                /* eliminate pairs */
                                                if (dup[p->sid] == -1) {
                                                    dup[p->sid] = p->sid;
                                                }
                                                else {
                                                    continue;
                                                }


                                                for (q = stNodes[n].lset[t]; q != NULL; q = q->next) {
                                                    f1 = p->sid;
                                                    f2 = q->sid;

                                                    if (f1 == f2) {
                                                        continue;
                                                    }

                                                    nPairs++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }


                /* merge the lsets of subtree */
                for (m = 0; m < SIGMA; m++) {
                    p = NULL;
                    for (j = sIndex + 1; j <= eIndex; j++) {
                        if ((q = stNodes[j].lset[m])) {

                            /* empty the subtree's ptrs array */
                            stNodes[j].lset[m] = NULL;
                            if (p == NULL) {
                                p = q;
                                stNodes[sIndex].lset[m] = q;
                            }
                            else {
                                p->next = q;
                            }

                            /* walk to the end */
                            while (p->next) {
                                p = p->next;
                            }
                        }
                    }
                }
            }
        }
        else {
            /* stnodes are sorted, so later part
             * will not satisfy EM cutoff */
            break;
        }
    }

    return nPairs;
}

