/*
 * $Rev: 77 $ 
 * $Date: 2011-05-11 16:03:49 -0700 (Wed, 11 May 2011) $ 
 * $Author: andy.cj.wu@gmail.com $
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * ----------------------------------------------------------------
 *
 */

#ifndef DYNAMIC_H_
#define DYNAMIC_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include "type.h"
#include "lib.h"
#include "elib.h"

#define NROW 2

#define CROW(i) ((i)%2)
#define PROW(i) ((i-1)%2)

#define BLOSUM62(map, ch1, ch2) \
        blosum62[map[ch1-'A']][map[ch2-'A']]

#define NEG_ADD(x, y) \
        (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y))

#define MAX(x, y) \
        (((x)>(y))? (x) : (y))

enum{OPEN = -10, GAP = -1};

#define MIN_VAL (-20000000)



/**
 * CELL for dynamic alignment result
 */
typedef struct cell{
    int score;           /* alignment score */
    int ndig;            /* #(matches) */
    int alen;            /* alignment length */
}CELL;

/**
 * parameters packed in struct
 */
typedef struct param{
    int AOL;             /* AlignOverLongerSeq */
    int SIM;             /* MatchSimilarity */
    int OS;              /* OptimalScoreOverSelfScore */
}PARAM;


void initMap(int nAA);
int selfScore(char *s, int ns);

CELL **allocTBL(int nrow, int ncol);
int **allocINT(int nrow, int ncol);
void freeTBL(CELL **tbl, int nrow);
int **allocINT(int nrow, int ncol);
void freeINT(int **tbl, int nrow);

void ffineGapAlign(char *s1, int s1Len, char *s2, int s2Len, CELL *result, CELL **tbl, int **del, int **ins);
void affineGapAlign(char *s1, int s1Len, char *s2, int s2Len, CELL *result, CELL **tbl, int **del, int **ins);
void printRow(CELL **tbl, int i, int ncol);

int isEdge(CELL *result, char *s1, int s1Len, char *s2, int s2Len, PARAM *param);

#endif /* DYNAMIC_H_ */


/*-------------------------------------------------------------------------------*
 *                             blosum62 matrix
 *                           -------------------    
       A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
    R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
    N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
    D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
    C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
    Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
    E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
    G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
    H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
    I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
    L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
    K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
    M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
    F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
    P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
    S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
    T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
    W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
    Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
    V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
    B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
    Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
    X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
    * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
*--------------------------------------------------------------------------------*/
