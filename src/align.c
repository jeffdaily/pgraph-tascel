/**
 * @file align.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "align.h"

#define CROW(i) ((i)%2)
#define PROW(i) ((i-1)%2)
#define SCORE(sub, map, first, ch1, ch2) \
    (sub)[(map)[(ch1)-(first)]][(map)[(ch2)-(first)]]
#define NEG_ADD(x, y) (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y))
#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))


pgraph_cell_t **pgraph_malloc_cell_table(size_t nrow, size_t ncol)
{
    size_t i = 0;
    pgraph_cell_t **table = NULL;

    table = (pgraph_cell_t **)malloc(nrow * sizeof(pgraph_cell_t *));
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        table[i] = malloc((ncol + 1) * sizeof(pgraph_cell_t));
    }

    return table;
}


int **pgraph_malloc_int_table(size_t nrow, size_t ncol)
{
    size_t i = 0;
    int **table = NULL;

    table = (int**)malloc(nrow * sizeof(int *));
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        table[i] = malloc((ncol + 1) * (sizeof * table[i]));
    }

    return table;
}


void pgraph_free_cell_table(pgraph_cell_t **table, size_t nrow)
{
    size_t i;
    for (i = 0; i < nrow; i++) {
        free(table[i]);
    }
    free(table);
}


void pgraph_free_int_table(int **table, size_t nrow)
{
    size_t i;
    for (i = 0; i < nrow; i++) {
        free(table[i]);
    }
    free(table);
}


int pgraph_self_score(
        const char * const restrict seq, size_t len,
        const int ** const restrict sub,
        const int * const restrict map, char first)
{
    size_t i = 0;
    int score = 0;

    for (i = 0; i < len; i++) {
        score += SCORE(sub, map, first, seq[i], seq[i]);
    }

    return score;
}


void pgraph_affine_gap_align(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        pgraph_cell_t * const restrict result,
        pgraph_cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        const int ** const restrict sub,
        const int * const restrict map, char first)
{
    size_t i;                                   /* first sequence char index */
    size_t j;                                   /* second sequence char index */
    pgraph_cell_t lastCol = {INT_MIN, 0, 0};    /* max cell in last col */
    pgraph_cell_t lastRow = {INT_MIN, 0, 0};    /* max cell in last row */
    pgraph_cell_t * restrict tI = NULL;         /* current DP table row */
    pgraph_cell_t * restrict pI = NULL;         /* previous DP table row */

    assert(s1Len > 0 && s2Len > 0);

    /* init first row of 3 tables */
    tbl[0][0].score = 0;
    tbl[0][0].matches = 0;
    tbl[0][0].length = 0;
    del[0][0] = 0;
    ins[0][0] = 0;
    tI = tbl[0];
    for (j = 1; j <= s2Len; j++) {
        tI[j].score = open + j * gap;
        tI[j].matches = 0;
        tI[j].length = 0;
        del[0][j] = INT_MIN;
        ins[0][j] = open + j * gap;
    }

    for (i = 1; i <= s1Len; ++i) {
        int dig = 0;        /* diagonal value in DP table */
        int up = 0;         /* upper value in DP table */
        int left = 0;       /* left value in DP table */
        int maxScore = 0;   /* max of dig, up, and left in DP table */
        size_t cr = 0;      /* current row index in DP table */
        size_t pr = 0;      /* previous row index in DP table */
        char ch1 = 0;       /* current character in first sequence */

        cr = CROW(i);
        pr = PROW(i);
        ch1 = s1[i - 1];
        tI = tbl[cr];
        pI = tbl[pr];

        /* init first column of 3 tables */
        tI[0].score = open + i * gap;
        tI[0].matches = 0;
        tI[0].length = 0;
        del[cr][0] = open + i * gap;
        ins[cr][0] = INT_MIN;

        for (j = 1; j <= s2Len; j++) {
            int tmp1 = 0;   /* temporary during DP calculation */
            int tmp2 = 0;   /* temporary during DP calculation */
            char ch2 = 0;   /* current character in second sequence */

            ch2 = s2[j - 1];

            /* overflow could happen, INT_MIN-1 = 2147483647 */
            tmp1 = pI[j].score + open + gap;
            tmp2 = NEG_ADD(del[pr][j], gap);
            up = MAX(tmp1, tmp2);
            del[cr][j] = up;
            tmp1 = tI[j - 1].score + open + gap;
            tmp2 = NEG_ADD(ins[cr][j - 1], gap);
            left = MAX(tmp1, tmp2);
            ins[cr][j] = left;
            maxScore = MAX(up, left);
            dig = pI[j - 1].score + SCORE(sub, map, first, ch1, ch2);
            maxScore = MAX(maxScore, dig);
            tI[j].score = maxScore;

            if (maxScore == dig) {
                tI[j].matches = pI[j - 1].matches + ((ch1 == ch2) ? 1 : 0);
                tI[j].length = pI[j - 1].length + 1;
            }
            else if (maxScore == up) {
                tI[j].matches = pI[j].matches;
                tI[j].length = pI[j].length + 1;
            }
            else {
                tI[j].matches = tI[j - 1].matches;
                tI[j].length = tI[j - 1].length + 1;
            }

            /* track the maximum of last row */
            if (i == s1Len && tI[j].score > lastRow.score) {
                lastRow = tI[j];
            }
        } /* end of j loop */

        assert(j == (s2Len + 1));

        /* update the maximum of last column */
        if (tI[s2Len].score > lastCol.score) {
            lastCol = tI[s2Len];
        }
    } /* end of i loop */

    *result = (lastCol.score > lastRow.score) ? lastCol : lastRow;
}


int pgraph_is_edge(
        const pgraph_cell_t result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int * const restrict self_score,
        size_t * const restrict max_len,
        const int ** const restrict sub,
        const int * const restrict map, char first)
{
    int sscore = 0;
    int maxLen = 0;
    int nmatch = 0;

    /* DO NOT need to compute the sscore value every time, if the later check
     * failed at the first step, the sscore computation is wasted */

    assert(s1_len);
    assert(s2_len);

    if (s1_len > s2_len) {
        maxLen = s1_len;
        sscore = pgraph_self_score(s1, s1_len, sub, map, first);
    }
    else {
        maxLen = s2_len;
        sscore = pgraph_self_score(s2, s2_len, sub, map, first);
    }

    nmatch = result.matches;
    *self_score = sscore;
    *max_len = maxLen;

    /* order the condition in strict->loose way, performance perspective
     * comparison using integers, no overflow could happen */
    if (result.score <= 0) {
        return 0;
    }
    else if ((result.length * 100 >= AOL * maxLen)
             && (nmatch * 100 >= SIM * result.length)
             && (result.score * 100 >= OS * sscore)) {
        return 1;
    }
    else {
        return 0;
    }
}

