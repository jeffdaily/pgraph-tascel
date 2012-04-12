/**
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef DYNAMIC_H_
#define DYNAMIC_H_

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#define SIGMA 26 /**< size of alphabet */
#define NROW 2
#define CROW(i) ((i)%2)
#define PROW(i) ((i-1)%2)
#define BLOSUM62(map, ch1, ch2) \
        blosum62[map[ch1-'A']][map[ch2-'A']]
#define NEG_ADD(x, y) \
        (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y))
#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))
#define MIN_VAL (-20000000)

enum {OPEN = -10, GAP = -1};
enum {NO = 0, YES = 1};
enum {FALSE = 0, TRUE = 1};

/**
 * cell_t for dynamic alignment result
 */
typedef struct {
    int score;           /* alignment score */
    int ndig;            /* #(matches) */
    int alen;            /* alignment length */
} cell_t;

/**
 * parameters for is_edge test, packed in struct
 */
typedef struct {
    int AOL;             /* AlignOverLongerSeq */
    int SIM;             /* MatchSimilarity */
    int OS;              /* OptimalScoreOverSelfScore */
} is_edge_param_t;


void init_map(int nAA);
int self_score(const char *s, size_t ns);

cell_t **alloc_tbl(int nrow, int ncol);
int **alloc_int(int nrow, int ncol);
void free_tbl(cell_t **tbl, int nrow);
int **alloc_int(int nrow, int ncol);
void free_int(int **tbl, int nrow);

void affine_gap_align(
        const char *s1, size_t s1Len,
        const char *s2, size_t s2Len,
        cell_t *result, cell_t **tbl, int **del, int **ins);

void print_row(cell_t **tbl, int i, int ncol);

int is_edge(const cell_t result,
        const char *s1, size_t s1Len,
        const char *s2, size_t s2Len,
        const is_edge_param_t param);

#ifdef __cplusplus
}
#endif 

#endif /* DYNAMIC_H_ */

/*----------------------------------------------------------------------------
 *                             blosum62 matrix
 * ---------------------------------------------------------------------------
 *     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
 *  A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
 *  R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
 *  N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
 *  D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 *  C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
 *  Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
 *  E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 *  G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
 *  H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
 *  I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
 *  L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
 *  K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
 *  M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
 *  F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
 *  P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 *  S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 *  T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
 *  W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
 *  Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 *  V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
 *  B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
 *  Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 *  X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
 *  * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
 * ---------------------------------------------------------------------------
 */
