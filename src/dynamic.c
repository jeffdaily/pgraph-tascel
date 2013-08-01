/**
 * @file dynamic.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "blosum/blosum40.h"
#include "blosum/blosum45.h"
#include "blosum/blosum62.h"
#include "blosum/blosum75.h"
#include "blosum/blosum80.h"
#include "blosum/blosum90.h"
#include "constants.h"
#include "dynamic.h"

#define CROW(i) ((i)%2)
#define PROW(i) ((i-1)%2)
#define BLOSUM(ch1, ch2) blosum[map[(ch1)-'A']][map[(ch2)-'A']]
#define NEG_ADD(x, y) (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y))
#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))
#define MIN_VAL (-20000000)

/* J stands for '*', 'O' and 'U' are dummy ones to make 26 */
static char AA[] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
                    'B', 'Z', 'X', 'J', 'O', 'U'
                   };

/** mapping from protein character to blosum index */
static int *map = NULL;

/** points to selected blosum table (default blosum62) */
static int (*blosum)[24] = blosum62;


/**
 * Initialize the protein index mapping used internally by the blosum matrix.
 */
static void init_blosum_map()
{
    if (NULL == map) {
        size_t i = 0;

        map = malloc(SIGMA * sizeof(int));
        for (i = 0; i < SIGMA; ++i) {
            map[AA[i] - 'A'] = (int)i;
        }
    }
}


int pg_match_blosum(char one, char two)
{
    return BLOSUM(one, two);
}


void pg_select_blosum(int number)
{
    init_blosum_map();

    switch (number) {
        case 40:
            blosum = blosum40;
            break;
        case 45:
            blosum = blosum45;
            break;
        case 62:
            blosum = blosum62;
            break;
        case 75:
            blosum = blosum75;
            break;
        case 80:
            blosum = blosum80;
            break;
        case 90:
            blosum = blosum90;
            break;
        default:
            fprintf(stderr, "invalid blosum number (%d)\n", number);
            assert(0);
    }
}


int pg_self_score(const sequence_t *s, match_t callback)
{
    size_t i = 0;
    int score = 0;

    for (i = 0; i < s->size; i++) {
        score += callback(s->str[i], s->str[i]);
    }

    return score;
}


int pg_self_score_blosum(const sequence_t *s)
{
    return pg_self_score_blosum2(s->str, s->size);
}


int pg_self_score_blosum2(const char *s, size_t len)
{
    size_t i = 0;
    int score = 0;

    init_blosum_map();

    for (i = 0; i < len; i++) {
        score += BLOSUM(s[i], s[i]);
    }

    return score;
}


cell_t **pg_alloc_tbl(int nrow, int ncol)
{
    int i;
    cell_t **tbl = NULL;

    assert(2 == nrow);

    tbl = (cell_t **)malloc(nrow * sizeof(cell_t *));
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        tbl[i] = malloc((ncol + 1) * sizeof(cell_t));
    }

    return tbl;
}


int **pg_alloc_int(int nrow, int ncol)
{
    int **tbl = NULL;
    int i;

    tbl = malloc(nrow * (sizeof * tbl));
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        tbl[i] = malloc((ncol + 1) * (sizeof * tbl[i]));
    }

    return tbl;
}


void pg_free_tbl(cell_t **tbl, int nrow)
{
    int i;
    for (i = 0; i < nrow; i++) {
        free(tbl[i]);
    }
    free(tbl);
}


void pg_free_int(int **tbl, int nrow)
{
    int i;
    for (i = 0; i < nrow; i++) {
        free(tbl[i]);
    }
    free(tbl);
}


void pg_affine_gap_align(
        const sequence_t *_s1, const sequence_t *_s2,
        cell_t *result, cell_t **tbl, int **del, int **ins,
        int open, int gap, match_t callback)
{
    size_t i;                           /* first sequence char index */
    size_t j;                           /* second sequence char index */
    const char *s1 = _s1->str;          /* char string of first sequence */
    size_t s1Len = _s1->size;           /* length of second sequence */
    const char *s2 = _s2->str;          /* char string of second sequence */
    size_t s2Len = _s2->size;           /* length of second sequence */
    cell_t lastCol = {INT_MIN, 0, 0};   /* max cell in last col */
    cell_t lastRow = {INT_MIN, 0, 0};   /* max cell in last row */
    cell_t *tI = NULL;                  /* current DP table row */
    cell_t *pI = NULL;                  /* previous DP table row */

    assert(s1Len > 0 && s2Len > 0);

    init_blosum_map();

    /* init first row of 3 tables */
    tbl[0][0].score = 0;
    tbl[0][0].ndig = 0;
    tbl[0][0].alen = 0;
    del[0][0] = 0;
    ins[0][0] = 0;
    tI = tbl[0];
    for (j = 1; j <= s2Len; j++) {
        tI[j].score = open + j * gap;
        tI[j].ndig = 0;
        tI[j].alen = 0;
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
        tI[0].ndig = 0;
        tI[0].alen = 0;
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
            dig = pI[j - 1].score + callback(ch1, ch2);
            maxScore = MAX(maxScore, dig);
            tI[j].score = maxScore;

            if (maxScore == dig) {
                tI[j].ndig = pI[j - 1].ndig + ((ch1 == ch2) ? 1 : 0);
                tI[j].alen = pI[j - 1].alen + 1;
            }
            else if (maxScore == up) {
                tI[j].ndig = pI[j].ndig;
                tI[j].alen = pI[j].alen + 1;
            }
            else {
                tI[j].ndig = tI[j - 1].ndig;
                tI[j].alen = tI[j - 1].alen + 1;
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


void pg_affine_gap_align_blosum2(
        const char *s1, size_t s1Len,
        const char *s2, size_t s2Len,
        cell_t *result, cell_t **tbl, int **del, int **ins,
        int open, int gap)
{
    size_t i;                           /* first sequence char index */
    size_t j;                           /* second sequence char index */
    cell_t lastCol = {INT_MIN, 0, 0};   /* max cell in last col */
    cell_t lastRow = {INT_MIN, 0, 0};   /* max cell in last row */
    cell_t *tI = NULL;                  /* current DP table row */
    cell_t *pI = NULL;                  /* previous DP table row */

    assert(s1Len > 0 && s2Len > 0);

    init_blosum_map();

    /* init first row of 3 tables */
    tbl[0][0].score = 0;
    tbl[0][0].ndig = 0;
    tbl[0][0].alen = 0;
    del[0][0] = 0;
    ins[0][0] = 0;
    tI = tbl[0];
    for (j = 1; j <= s2Len; j++) {
        tI[j].score = open + j * gap;
        tI[j].ndig = 0;
        tI[j].alen = 0;
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
        tI[0].ndig = 0;
        tI[0].alen = 0;
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
            dig = pI[j - 1].score + BLOSUM(ch1, ch2);
            maxScore = MAX(maxScore, dig);
            tI[j].score = maxScore;

            if (maxScore == dig) {
                tI[j].ndig = pI[j - 1].ndig + ((ch1 == ch2) ? 1 : 0);
                tI[j].alen = pI[j - 1].alen + 1;
            }
            else if (maxScore == up) {
                tI[j].ndig = pI[j].ndig;
                tI[j].alen = pI[j].alen + 1;
            }
            else {
                tI[j].ndig = tI[j - 1].ndig;
                tI[j].alen = tI[j - 1].alen + 1;
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


void pg_affine_gap_align_blosum(
        const sequence_t *s1, const sequence_t *s2,
        cell_t *result, cell_t **tbl, int **del, int **ins,
        int open, int gap)
{
    pg_affine_gap_align_blosum2(s1->str, s1->size, s2->str, s2->size, result, tbl, del, ins, open, gap);
}


void pg_print_row(cell_t **tbl, int i, int ncol)
{
    int j;
    printf("[");
    for (j = 0; j <= ncol; j++) {
        printf("%d, %d\t", tbl[i][j].score, tbl[i][j].ndig);
    }
    printf("]\n");
}


int pg_is_edge(const cell_t result,
               const sequence_t *s1, const sequence_t *s2,
               const param_t param, int *_sscore, size_t *_maxLen,
               match_t callback)
{
    int sscore = 0;
    int maxLen = 0;
    int nmatch = 0;

    /* DO NOT need to compute the sscore value every time, if the later check
     * failed at the first step, the sscore computation is wasted */

    assert(s1->size);
    assert(s2->size);

    if (s1->size > s2->size) {
        maxLen = s1->size;
        sscore = pg_self_score(s1, callback);
    }
    else {
        maxLen = s2->size;
        sscore = pg_self_score(s2, callback);
    }

    nmatch = result.ndig;
    *_sscore = sscore;
    *_maxLen = maxLen;

    /* order the condition in strict->loose way, performance perspective
     * comparison using integers, no overflow could happen */
    if (result.score <= 0) {
        return FALSE;
    }
    else if ((result.alen * 100 >= param.AOL * maxLen)
             && (nmatch * 100 >= param.SIM * result.alen)
             && (result.score * 100 >= param.OS * sscore)) {
        return TRUE;
    }
    else {
        return FALSE;
    }
}


int pg_is_edge_blosum(const cell_t result,
               const sequence_t *s1, const sequence_t *s2,
               const param_t param, int *_sscore, size_t *_maxLen)
{
    int sscore = 0;
    int maxLen = 0;
    int nmatch = 0;

    /* DO NOT need to compute the sscore value every time, if the later check
     * failed at the first step, the sscore computation is wasted */

    assert(s1->size);
    assert(s2->size);

    if (s1->size > s2->size) {
        maxLen = s1->size;
        sscore = pg_self_score_blosum(s1);
    }
    else {
        maxLen = s2->size;
        sscore = pg_self_score_blosum(s2);
    }

    nmatch = result.ndig;
    *_sscore = sscore;
    *_maxLen = maxLen;

    /* order the condition in strict->loose way, performance perspective
     * comparison using integers, no overflow could happen */
    if (result.score <= 0) {
        return FALSE;
    }
    else if ((result.alen * 100 >= param.AOL * maxLen)
             && (nmatch * 100 >= param.SIM * result.alen)
             && (result.score * 100 >= param.OS * sscore)) {
        return TRUE;
    }
    else {
        return FALSE;
    }
}

