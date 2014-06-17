/**
 * @file alignment.cpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "alignment.hpp"
#include "blosum/blosum40.h"
#include "blosum/blosum45.h"
#include "blosum/blosum50.h"
#include "blosum/blosum62.h"
#include "blosum/blosum75.h"
#include "blosum/blosum80.h"
#include "blosum/blosum90.h"
#if HAVE_EMMINTRIN_H
#include "ssw.h"

extern "C" {
#include "global_sse2.h"
#include "defs.h"
#include "param.h"
#include "dropgsw2.h"
}
#endif

/** mapping from BLOSUM alphabet to BLOSUM index; use as BLOSUM[map[ch-'A']] */
static const int MAP_BLOSUM[] = {
    0,  20,   4,   3,   6,  13,   7,   8,   9,  23,
    11,  10,  12,   2,  23,  14,   5,   1,  15,  16,
    22,  19,  17,  22,  18,  21 };

/* This table is used to transform amino acid letters into numbers. */
static const int MAP_BLOSUM_[128] = {
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 22, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 22, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};


#define CROW(i) ((i)%2)
#define PROW(i) ((i-1)%2)
//#define BLOSUM(ch1, ch2) (blosum[MAP_BLOSUM[(ch1)-'A']][MAP_BLOSUM[(ch2)-'A']])
#define BLOSUM(ch1, ch2) (blosum[MAP_BLOSUM_[(ch1)]][MAP_BLOSUM_[(ch2)]])
#define SCORE(sub, map, first, ch1, ch2) \
    (sub)[(map)[(ch1)-(first)]][(map)[(ch2)-(first)]]
#define NEG_ADD(x, y) (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y))
#define MAX(x, y) (((x)>(y))? (x) : (y))
//#define NEG_INF (INT_MIN/2)
static const int NEG_INF = SHRT_MIN / 2;


using std::cerr;
using std::endl;


namespace pgraph {


/** points to selected blosum table (default blosum62) */
static const int (* restrict blosum)[24] = blosum62;
static const int * const restrict * restrict blosum_ = blosum62_;
static const int8_t * restrict blosum__ = blosum62__;


cell_t **allocate_cell_table(size_t nrow, size_t ncol)
{
    return allocate_table<cell_t>(nrow, ncol);
}


tbl_t **allocate_tbl_table(size_t nrow, size_t ncol)
{
    return allocate_table<tbl_t>(nrow, ncol);
}


int **allocate_int_table(size_t nrow, size_t ncol)
{
    return allocate_table<int>(nrow, ncol);
}


void free_cell_table(cell_t **table, size_t nrow)
{
    free_table<cell_t>(table, nrow);
}


void free_int_table(int **table, size_t nrow)
{
    free_table<int>(table, nrow);
}


int self_score(
        const char * const restrict seq, size_t len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first)
{
    size_t i = 0;
    int score = 0;

    for (i = 0; i < len; i++) {
        score += SCORE(sub, map, first, seq[i], seq[i]);
    }

    return score;
}


int self_score(const char * const restrict seq, size_t len, match_t callback)
{
    size_t i = 0;
    int score = 0;

    for (i = 0; i < len; i++) {
        score += callback(seq[i], seq[i]);
    }

    return score;
}


int self_score(const char * const restrict /*seq*/, size_t len, int match)
{
    return match * len;
}


int self_score_blosum(const char * const restrict seq, size_t len)
{
    size_t i = 0;
    int score = 0;

    for (i = 0; i < len; i++) {
        score += BLOSUM(seq[i], seq[i]);
    }

    return score;
}


void select_blosum(int number)
{
    switch (number) {
        case 40:
            blosum = blosum40;
            blosum_ = blosum40_;
            blosum__ = blosum40__;
            break;
        case 45:
            blosum = blosum45;
            blosum_ = blosum45_;
            blosum__ = blosum45__;
            break;
        case 62:
            blosum = blosum62;
            blosum_ = blosum62_;
            blosum__ = blosum62__;
            break;
        case 75:
            blosum = blosum75;
            blosum_ = blosum75_;
            blosum__ = blosum75__;
            break;
        case 80:
            blosum = blosum80;
            blosum_ = blosum80_;
            blosum__ = blosum80__;
            break;
        case 90:
            blosum = blosum90;
            blosum_ = blosum90_;
            blosum__ = blosum90__;
            break;
        default:
            cerr << "invalid blosum number (" << number << ")" << endl;
            assert(0);
    }
}

#if USE_SIMILARITIES

#define AFFINE_GAP_DECL                                                     \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    bool delete_tbl = (tbl_ == NULL);                                       \
    bool delete_del = (del_ == NULL);                                       \
    bool delete_ins = (ins_ == NULL);                                       \
    size_t longestLen = MAX(s1Len,s2Len);                                   \
    cell_t * const restrict * const restrict tbl =                          \
            delete_tbl ? allocate_cell_table(2U,longestLen) : tbl_;         \
    int * const restrict * const restrict del =                             \
            delete_del ? allocate_int_table(2U,longestLen) : del_;          \
    int * const restrict * const restrict ins =                             \
            delete_ins ? allocate_int_table(2U,longestLen) : ins_;          \
    cell_t maxCell = {NEG_INF, 0, 0, 0};   /* max cell in table */          \

#define AFFINE_GAP_DECL_SEMI                                                \
    cell_t lastCol = {NEG_INF, 0, 0, 0};   /* max cell in last col */       \
    cell_t lastRow = {NEG_INF, 0, 0, 0};   /* max cell in last row */       \

#define AFFINE_GAP_ASSERT       \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(tbl);                \
    assert(del);                \
    assert(ins);

#define AFFINE_GAP_INIT_CORNER                                              \
    tbl[0][0].score = 0;                                                    \
    tbl[0][0].matches = 0;                                                  \
    tbl[0][0].similarities = 0;                                             \
    tbl[0][0].length = 0;                                                   \
    del[0][0] = NEG_INF;                                                    \
    ins[0][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY                              \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = open + (j-1) * gap;                               \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].similarities = 0;                                         \
        tbl[0][j].length = 0;                                               \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY                           \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = 0;                                                \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].similarities = 0;                                         \
        tbl[0][j].length = 0;                                               \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_BODY1                                                    \
    for (i = 1; i <= s1Len; ++i) {                                          \
        size_t cr = CROW(i);                                                \
        size_t pr = PROW(i);                                                \
        char ch1 = s1[i - 1];                                               \
        const int * const restrict BlosumRow = blosum_[MAP_BLOSUM_[ch1]];   \
                                                                            \
        int Nscore        = tbl[pr][0].score;                               \
        int Nmatches      = tbl[pr][0].matches;                             \
        int Nsimilarities = tbl[pr][0].similarities;                        \
        int Nlength       = tbl[pr][0].length;                              \
                                                                            \
        int Wscore        = open + (i-1) * gap;                             \
        int Wmatches      = 0;                                              \
        int Wsimilarities = 0;                                              \
        int Wlength       = 0;

#define AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY                              \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = Wscore;                                          \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].similarities = Wsimilarities;                            \
        tbl[cr][0].length = Wlength;                                        \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY                           \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = 0;                                               \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].similarities = Wsimilarities;                            \
        tbl[cr][0].length = Wlength;                                        \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_BODY2(CALCULATE_SCORE)                                   \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int up = 0;                                                     \
            int left = 0;                                                   \
            int dig = 0;                                                    \
            char ch2 = s2[j - 1];                                           \
            int NWscore        = Nscore;                                    \
            int NWmatches      = Nmatches;                                  \
            int NWsimilarities = Nsimilarities;                             \
            int NWlength       = Nlength;                                   \
                                                                            \
            Nscore        = tbl[pr][j].score;                               \
            Nmatches      = tbl[pr][j].matches;                             \
            Nsimilarities = tbl[pr][j].similarities;                        \
            Nlength       = tbl[pr][j].length;                              \
                                                                            \
            tmp1 = CALCULATE_SCORE;                                         \
            up   = MAX(Nscore + open, del[pr][j]   + gap);                  \
            left = MAX(Wscore + open, ins[cr][j-1] + gap);                  \
            dig  = NWscore + tmp1;                                          \
                                                                            \
            del[cr][j] = up;                                                \
            ins[cr][j] = left;                                              \
                                                                            \
            if ((dig >= up) && (dig >= left)) {                             \
                Wscore = dig;                                               \
                Wmatches  = NWmatches + (ch1 == ch2);                       \
                Wsimilarities = NWsimilarities + (tmp1 > 0);                \
                Wlength  = NWlength + 1;                                    \
            } else if (up >= left) {                                        \
                Wscore = up;                                                \
                Wmatches  = Nmatches;                                       \
                Wsimilarities = Nsimilarities;                              \
                Wlength  = Nlength + 1;                                     \
            } else {                                                        \
                Wscore = left;                                              \
                Wmatches  = Wmatches;                                       \
                Wsimilarities  = Wsimilarities;                             \
                Wlength  = Wlength + 1;                                     \
            }                                                               \
                                                                            \
            tbl[cr][j].score = Wscore;                                      \
            tbl[cr][j].matches  = Wmatches;                                 \
            tbl[cr][j].similarities  = Wsimilarities;                       \
            tbl[cr][j].length  = Wlength;

#define AFFINE_GAP_BODY2_LOCAL                                              \
            if (tbl[cr][j].score > maxCell.score) {                         \
                maxCell = tbl[cr][j];                                       \
            }                                                               \
            if (tbl[cr][j].score <= 0) {                                    \
                tbl[cr][j].score = Wscore = 0;                              \
                tbl[cr][j].matches = Wmatches = 0;                          \
                tbl[cr][j].similarities = Wsimilarities = 0;                \
                tbl[cr][j].length = Wlength = 0;                            \
            }

#define AFFINE_GAP_END_GLOBAL                                               \
        } /* end of j loop */                                               \
    } /* end of i loop */                                                   \
    maxCell = tbl[CROW(s1Len)][s2Len];

#define AFFINE_GAP_END_SEMI                                                 \
            /* track the maximum of last row */                             \
            if (i == s1Len && tbl[cr][j].score > lastRow.score) {           \
                lastRow = tbl[cr][j];                                       \
            }                                                               \
        } /* end of j loop */                                               \
        /* update the maximum of last column */                             \
        if (tbl[cr][s2Len].score > lastCol.score) {                         \
            lastCol = tbl[cr][s2Len];                                       \
        }                                                                   \
    } /* end of i loop */                                                   \
    maxCell = (lastCol.score > lastRow.score) ? lastCol : lastRow;

#define AFFINE_GAP_END_LOCAL                                                \
        } /* end of j loop */                                               \
    } /* end of i loop */

#define AFFINE_GAP_FREE                                                     \
    if (delete_tbl) free_table(tbl,2U);                                     \
    if (delete_del) free_table(del,2U);                                     \
    if (delete_ins) free_table(ins,2U);

#define AFFINE_GAP_RETURN                                                   \
    return maxCell;

#else /* USE_SIMILARITIES */

#define AFFINE_GAP_DECL                                                     \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    bool delete_tbl = (tbl_ == NULL);                                       \
    bool delete_del = (del_ == NULL);                                       \
    bool delete_ins = (ins_ == NULL);                                       \
    size_t longestLen = MAX(s1Len,s2Len);                                   \
    cell_t * const restrict * const restrict tbl =                          \
            delete_tbl ? allocate_cell_table(2U,longestLen) : tbl_;         \
    int * const restrict * const restrict del =                             \
            delete_del ? allocate_int_table(2U,longestLen) : del_;          \
    int * const restrict * const restrict ins =                             \
            delete_ins ? allocate_int_table(2U,longestLen) : ins_;          \
    cell_t maxCell = {NEG_INF, 0, 0};   /* max cell in table */             \

#define AFFINE_GAP_DECL_SEMI                                                \
    cell_t lastCol = {NEG_INF, 0, 0};   /* max cell in last col */          \
    cell_t lastRow = {NEG_INF, 0, 0};   /* max cell in last row */          \

#define AFFINE_GAP_ASSERT       \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(tbl);                \
    assert(del);                \
    assert(ins);

#define AFFINE_GAP_INIT_CORNER                                              \
    tbl[0][0].score = 0;                                                    \
    tbl[0][0].matches = 0;                                                  \
    tbl[0][0].length = 0;                                                   \
    del[0][0] = NEG_INF;                                                    \
    ins[0][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY                              \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = open + (j-1) * gap;                               \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].length = 0;                                               \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY                           \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = 0;                                                \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].length = 0;                                               \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_BODY1                                                    \
    for (i = 1; i <= s1Len; ++i) {                                          \
        size_t cr = CROW(i);                                                \
        size_t pr = PROW(i);                                                \
        char ch1 = s1[i - 1];                                               \
        const int * const restrict BlosumRow = blosum_[MAP_BLOSUM_[ch1]];   \
                                                                            \
        int Nscore        = tbl[pr][0].score;                               \
        int Nmatches      = tbl[pr][0].matches;                             \
        int Nlength       = tbl[pr][0].length;                              \
                                                                            \
        int Wscore        = open + (i-1) * gap;                             \
        int Wmatches      = 0;                                              \
        int Wlength       = 0;

#define AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY                              \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = Wscore;                                          \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].length = Wlength;                                        \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY                           \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = 0;                                               \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].length = Wlength;                                        \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_BODY2(CALCULATE_SCORE)                                   \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int up = 0;                                                     \
            int left = 0;                                                   \
            int dig = 0;                                                    \
            char ch2 = s2[j - 1];                                           \
            int NWscore        = Nscore;                                    \
            int NWmatches      = Nmatches;                                  \
            int NWlength       = Nlength;                                   \
                                                                            \
            Nscore        = tbl[pr][j].score;                               \
            Nmatches      = tbl[pr][j].matches;                             \
            Nlength       = tbl[pr][j].length;                              \
                                                                            \
            tmp1 = CALCULATE_SCORE;                                         \
            up   = MAX(Nscore + open, del[pr][j]   + gap);                  \
            left = MAX(Wscore + open, ins[cr][j-1] + gap);                  \
            dig  = NWscore + tmp1;                                          \
                                                                            \
            del[cr][j] = up;                                                \
            ins[cr][j] = left;                                              \
                                                                            \
            if ((dig >= up) && (dig >= left)) {                             \
                Wscore = dig;                                               \
                Wmatches  = NWmatches + (ch1 == ch2);                       \
                Wlength  = NWlength + 1;                                    \
            } else if (up >= left) {                                        \
                Wscore = up;                                                \
                Wmatches  = Nmatches;                                       \
                Wlength  = Nlength + 1;                                     \
            } else {                                                        \
                Wscore = left;                                              \
                Wmatches  = Wmatches;                                       \
                Wlength  = Wlength + 1;                                     \
            }                                                               \
                                                                            \
            tbl[cr][j].score = Wscore;                                      \
            tbl[cr][j].matches  = Wmatches;                                 \
            tbl[cr][j].length  = Wlength;

#define AFFINE_GAP_BODY2_LOCAL                                              \
            if (tbl[cr][j].score > maxCell.score) {                         \
                maxCell = tbl[cr][j];                                       \
            }                                                               \
            if (tbl[cr][j].score <= 0) {                                    \
                tbl[cr][j].score = Wscore = 0;                              \
                tbl[cr][j].matches = Wmatches = 0;                          \
                tbl[cr][j].length = Wlength = 0;                            \
            }

#define AFFINE_GAP_END_GLOBAL                                               \
        } /* end of j loop */                                               \
    } /* end of i loop */                                                   \
    maxCell = tbl[CROW(s1Len)][s2Len];

#define AFFINE_GAP_END_SEMI                                                 \
            /* track the maximum of last row */                             \
            if (i == s1Len && tbl[cr][j].score > lastRow.score) {           \
                lastRow = tbl[cr][j];                                       \
            }                                                               \
        } /* end of j loop */                                               \
        /* update the maximum of last column */                             \
        if (tbl[cr][s2Len].score > lastCol.score) {                         \
            lastCol = tbl[cr][s2Len];                                       \
        }                                                                   \
    } /* end of i loop */                                                   \
    maxCell = (lastCol.score > lastRow.score) ? lastCol : lastRow;

#define AFFINE_GAP_END_LOCAL                                                \
        } /* end of j loop */                                               \
    } /* end of i loop */

#define AFFINE_GAP_FREE                                                     \
    if (delete_tbl) free_table(tbl,2U);                                     \
    if (delete_del) free_table(del,2U);                                     \
    if (delete_ins) free_table(ins,2U);

#define AFFINE_GAP_RETURN                                                   \
    return maxCell;

#endif /* USE_SIMILARITIES */



cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        cell_t * const restrict * const restrict tbl_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


#undef AFFINE_GAP_DECL
#undef AFFINE_GAP_DECL_SEMI
#undef AFFINE_GAP_ASSERT
#undef AFFINE_GAP_INIT_CORNER
#undef AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
#undef AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
#undef AFFINE_GAP_BODY1
#undef AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
#undef AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
#undef AFFINE_GAP_BODY2
#undef AFFINE_GAP_BODY2_LOCAL
#undef AFFINE_GAP_END_GLOBAL
#undef AFFINE_GAP_END_SEMI
#undef AFFINE_GAP_END_LOCAL
#undef AFFINE_GAP_FREE
#undef AFFINE_GAP_RETURN


#if USE_SIMILARITIES

#define AFFINE_GAP_DECL                                                     \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    bool delete_tbl = (tbl_ == NULL);                                       \
    size_t longestLen = MAX(s1Len,s2Len);                                   \
    tbl_t * const restrict * const restrict tbl =                           \
            delete_tbl ? allocate_tbl_table(2U,longestLen) : tbl_;          \
    tbl_t maxCell = {NEG_INF, 0, 0, 0, 0, 0};   /* max cell in table */

#define AFFINE_GAP_DECL_SEMI                                                \
    tbl_t lastCol = {NEG_INF, 0, 0, 0, 0, 0};   /* max cell in last col */ \
    tbl_t lastRow = {NEG_INF, 0, 0, 0, 0, 0};   /* max cell in last row */

#define AFFINE_GAP_ASSERT       \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(tbl);

#define AFFINE_GAP_INIT_CORNER                                              \
    tbl[0][0].score = 0;                                                    \
    tbl[0][0].matches = 0;                                                  \
    tbl[0][0].similarities = 0;                                             \
    tbl[0][0].length = 0;                                                   \
    tbl[0][0].del = NEG_INF;                                                    \
    tbl[0][0].ins = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY                              \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = open + (j-1) * gap;                               \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].similarities = 0;                                         \
        tbl[0][j].length = 0;                                               \
        tbl[0][j].del = NEG_INF;                                                \
        tbl[0][j].ins = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY                           \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = 0;                                                \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].similarities = 0;                                         \
        tbl[0][j].length = 0;                                               \
        tbl[0][j].del = NEG_INF;                                                \
        tbl[0][j].ins = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_BODY1                                                    \
    for (i = 1; i <= s1Len; ++i) {                                          \
        size_t cr = CROW(i);                                                \
        size_t pr = PROW(i);                                                \
        char ch1 = s1[i - 1];                                               \
        const int * const restrict BlosumRow = blosum_[MAP_BLOSUM_[ch1]];   \
                                                                            \
        int Nscore        = tbl[pr][0].score;                               \
        int Nmatches      = tbl[pr][0].matches;                             \
        int Nsimilarities = tbl[pr][0].similarities;                        \
        int Nlength       = tbl[pr][0].length;                              \
                                                                            \
        int Wscore        = open + (i-1) * gap;                             \
        int Wmatches      = 0;                                              \
        int Wsimilarities = 0;                                              \
        int Wlength       = 0;

#define AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY                              \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = Wscore;                                          \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].similarities = Wsimilarities;                            \
        tbl[cr][0].length = Wlength;                                        \
        tbl[cr][0].del = open + (i-1) * gap;                                    \
        tbl[cr][0].ins = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY                           \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = 0;                                               \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].similarities = Wsimilarities;                            \
        tbl[cr][0].length = Wlength;                                        \
        tbl[cr][0].del = open + (i-1) * gap;                                    \
        tbl[cr][0].ins = NEG_INF;

#define AFFINE_GAP_BODY2(CALCULATE_SCORE)                                   \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int up = 0;                                                     \
            int left = 0;                                                   \
            int dig = 0;                                                    \
            char ch2 = s2[j - 1];                                           \
            int NWscore        = Nscore;                                    \
            int NWmatches      = Nmatches;                                  \
            int NWsimilarities = Nsimilarities;                             \
            int NWlength       = Nlength;                                   \
                                                                            \
            Nscore        = tbl[pr][j].score;                               \
            Nmatches      = tbl[pr][j].matches;                             \
            Nsimilarities = tbl[pr][j].similarities;                        \
            Nlength       = tbl[pr][j].length;                              \
                                                                            \
            tmp1 = CALCULATE_SCORE;                                         \
            up   = MAX(Nscore + open, tbl[pr][j].del   + gap);                  \
            left = MAX(Wscore + open, tbl[cr][j-1].ins + gap);                  \
            dig  = NWscore + tmp1;                                          \
                                                                            \
            tbl[cr][j].del = up;                                                \
            tbl[cr][j].ins = left;                                              \
                                                                            \
            if ((dig >= up) && (dig >= left)) {                             \
                Wscore = dig;                                               \
                Wmatches  = NWmatches + (ch1 == ch2);                       \
                Wsimilarities = NWsimilarities + (tmp1 > 0);                \
                Wlength  = NWlength + 1;                                    \
            } else if (up >= left) {                                        \
                Wscore = up;                                                \
                Wmatches  = Nmatches;                                       \
                Wsimilarities = Nsimilarities;                              \
                Wlength  = Nlength + 1;                                     \
            } else {                                                        \
                Wscore = left;                                              \
                Wmatches  = Wmatches;                                       \
                Wsimilarities  = Wsimilarities;                             \
                Wlength  = Wlength + 1;                                     \
            }                                                               \
                                                                            \
            tbl[cr][j].score = Wscore;                                      \
            tbl[cr][j].matches  = Wmatches;                                 \
            tbl[cr][j].similarities  = Wsimilarities;                       \
            tbl[cr][j].length  = Wlength;

#define AFFINE_GAP_BODY2_LOCAL                                              \
            if (tbl[cr][j].score > maxCell.score) {                         \
                maxCell = tbl[cr][j];                                       \
            }                                                               \
            if (tbl[cr][j].score <= 0) {                                    \
                tbl[cr][j].score = Wscore = 0;                              \
                tbl[cr][j].matches = Wmatches = 0;                          \
                tbl[cr][j].similarities = Wsimilarities = 0;                \
                tbl[cr][j].length = Wlength = 0;                            \
            }

#define AFFINE_GAP_END_GLOBAL                                               \
        } /* end of j loop */                                               \
    } /* end of i loop */                                                   \
    maxCell = tbl[CROW(s1Len)][s2Len];

#define AFFINE_GAP_END_SEMI                                                 \
            /* track the maximum of last row */                             \
            if (i == s1Len && tbl[cr][j].score > lastRow.score) {           \
                lastRow = tbl[cr][j];                                       \
            }                                                               \
        } /* end of j loop */                                               \
        /* update the maximum of last column */                             \
        if (tbl[cr][s2Len].score > lastCol.score) {                         \
            lastCol = tbl[cr][s2Len];                                       \
        }                                                                   \
    } /* end of i loop */                                                   \
    maxCell = (lastCol.score > lastRow.score) ? lastCol : lastRow;

#define AFFINE_GAP_END_LOCAL                                                \
        } /* end of j loop */                                               \
    } /* end of i loop */

#define AFFINE_GAP_FREE                                                     \
    if (delete_tbl) free_table(tbl,2U);

#define AFFINE_GAP_RETURN                                                   \
    cell_t retval = {                                                       \
        maxCell.score,                                                      \
        maxCell.matches,                                                    \
        maxCell.similarities,                                               \
        maxCell.length                                                      \
    };                                                                      \
    return retval;

#else /* USE_SIMILARITIES */

#define AFFINE_GAP_DECL                                                     \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    bool delete_tbl = (tbl_ == NULL);                                       \
    size_t longestLen = MAX(s1Len,s2Len);                                   \
    tbl_t * const restrict * const restrict tbl =                           \
            delete_tbl ? allocate_tbl_table(2U,longestLen) : tbl_;          \
    tbl_t maxCell = {NEG_INF, 0, 0, 0, 0};   /* max cell in table */

#define AFFINE_GAP_DECL_SEMI                                                \
    tbl_t lastCol = {NEG_INF, 0, 0, 0, 0};   /* max cell in last col */     \
    tbl_t lastRow = {NEG_INF, 0, 0, 0, 0};   /* max cell in last row */

#define AFFINE_GAP_ASSERT       \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(tbl);

#define AFFINE_GAP_INIT_CORNER                                              \
    tbl[0][0].score = 0;                                                    \
    tbl[0][0].matches = 0;                                                  \
    tbl[0][0].length = 0;                                                   \
    tbl[0][0].del = NEG_INF;                                                    \
    tbl[0][0].ins = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY                              \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = open + (j-1) * gap;                               \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].length = 0;                                               \
        tbl[0][j].del = NEG_INF;                                                \
        tbl[0][j].ins = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY                           \
    for (j = 1; j <= s2Len; j++) {                                          \
        tbl[0][j].score = 0;                                                \
        tbl[0][j].matches = 0;                                              \
        tbl[0][j].length = 0;                                               \
        tbl[0][j].del = NEG_INF;                                                \
        tbl[0][j].ins = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_BODY1                                                    \
    for (i = 1; i <= s1Len; ++i) {                                          \
        size_t cr = CROW(i);                                                \
        size_t pr = PROW(i);                                                \
        char ch1 = s1[i - 1];                                               \
        const int * const restrict BlosumRow = blosum_[MAP_BLOSUM_[ch1]];   \
                                                                            \
        int Nscore        = tbl[pr][0].score;                               \
        int Nmatches      = tbl[pr][0].matches;                             \
        int Nlength       = tbl[pr][0].length;                              \
                                                                            \
        int Wscore        = open + (i-1) * gap;                             \
        int Wmatches      = 0;                                              \
        int Wlength       = 0;

#define AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY                              \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = Wscore;                                          \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].length = Wlength;                                        \
        tbl[cr][0].del = open + (i-1) * gap;                                    \
        tbl[cr][0].ins = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY                           \
        /* init first column of 3 tables */                                 \
        tbl[cr][0].score = 0;                                               \
        tbl[cr][0].matches = Wmatches;                                      \
        tbl[cr][0].length = Wlength;                                        \
        tbl[cr][0].del = open + (i-1) * gap;                                    \
        tbl[cr][0].ins = NEG_INF;

#define AFFINE_GAP_BODY2(CALCULATE_SCORE)                                   \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int up = 0;                                                     \
            int left = 0;                                                   \
            int dig = 0;                                                    \
            char ch2 = s2[j - 1];                                           \
            int NWscore        = Nscore;                                    \
            int NWmatches      = Nmatches;                                  \
            int NWlength       = Nlength;                                   \
                                                                            \
            Nscore        = tbl[pr][j].score;                               \
            Nmatches      = tbl[pr][j].matches;                             \
            Nlength       = tbl[pr][j].length;                              \
                                                                            \
            tmp1 = CALCULATE_SCORE;                                         \
            up   = MAX(Nscore + open, tbl[pr][j].del   + gap);                  \
            left = MAX(Wscore + open, tbl[cr][j-1].ins + gap);                  \
            dig  = NWscore + tmp1;                                          \
                                                                            \
            tbl[cr][j].del = up;                                                \
            tbl[cr][j].ins = left;                                              \
                                                                            \
            if ((dig >= up) && (dig >= left)) {                             \
                Wscore = dig;                                               \
                Wmatches  = NWmatches + (ch1 == ch2);                       \
                Wlength  = NWlength + 1;                                    \
            } else if (up >= left) {                                        \
                Wscore = up;                                                \
                Wmatches  = Nmatches;                                       \
                Wlength  = Nlength + 1;                                     \
            } else {                                                        \
                Wscore = left;                                              \
                Wmatches  = Wmatches;                                       \
                Wlength  = Wlength + 1;                                     \
            }                                                               \
                                                                            \
            tbl[cr][j].score = Wscore;                                      \
            tbl[cr][j].matches  = Wmatches;                                 \
            tbl[cr][j].length  = Wlength;

#define AFFINE_GAP_BODY2_LOCAL                                              \
            if (tbl[cr][j].score > maxCell.score) {                         \
                maxCell = tbl[cr][j];                                       \
            }                                                               \
            if (tbl[cr][j].score <= 0) {                                    \
                tbl[cr][j].score = Wscore = 0;                              \
                tbl[cr][j].matches = Wmatches = 0;                          \
                tbl[cr][j].length = Wlength = 0;                            \
            }

#define AFFINE_GAP_END_GLOBAL                                               \
        } /* end of j loop */                                               \
    } /* end of i loop */                                                   \
    maxCell = tbl[CROW(s1Len)][s2Len];

#define AFFINE_GAP_END_SEMI                                                 \
            /* track the maximum of last row */                             \
            if (i == s1Len && tbl[cr][j].score > lastRow.score) {           \
                lastRow = tbl[cr][j];                                       \
            }                                                               \
        } /* end of j loop */                                               \
        /* update the maximum of last column */                             \
        if (tbl[cr][s2Len].score > lastCol.score) {                         \
            lastCol = tbl[cr][s2Len];                                       \
        }                                                                   \
    } /* end of i loop */                                                   \
    maxCell = (lastCol.score > lastRow.score) ? lastCol : lastRow;

#define AFFINE_GAP_END_LOCAL                                                \
        } /* end of j loop */                                               \
    } /* end of i loop */

#define AFFINE_GAP_FREE                                                     \
    if (delete_tbl) free_table(tbl,2U);

#define AFFINE_GAP_RETURN                                                   \
    cell_t retval = {                                                       \
        maxCell.score,                                                      \
        maxCell.matches,                                                    \
        maxCell.length                                                      \
    };                                                                      \
    return retval;

#endif /* USE_SIMILARITIES */



cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        tbl_t * const restrict * const restrict tbl_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


#undef AFFINE_GAP_DECL
#undef AFFINE_GAP_DECL_SEMI
#undef AFFINE_GAP_ASSERT
#undef AFFINE_GAP_INIT_CORNER
#undef AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
#undef AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
#undef AFFINE_GAP_BODY1
#undef AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
#undef AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
#undef AFFINE_GAP_BODY2
#undef AFFINE_GAP_BODY2_LOCAL
#undef AFFINE_GAP_END_GLOBAL
#undef AFFINE_GAP_END_SEMI
#undef AFFINE_GAP_END_LOCAL
#undef AFFINE_GAP_FREE
#undef AFFINE_GAP_RETURN


#if USE_SIMILARITIES

#define AFFINE_GAP_DECL                                                     \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    bool delete_scr = (scr_ == NULL);                                       \
    bool delete_mat = (mat_ == NULL);                                       \
    bool delete_sim = (sim_ == NULL);                                       \
    bool delete_len = (len_ == NULL);                                       \
    bool delete_del = (del_ == NULL);                                       \
    bool delete_ins = (ins_ == NULL);                                       \
    size_t longestLen = MAX(s1Len,s2Len);                                   \
    int * const restrict * const restrict scr =                             \
            delete_scr ? allocate_int_table(2U,longestLen) : scr_;          \
    int * const restrict * const restrict mat =                             \
            delete_mat ? allocate_int_table(2U,longestLen) : mat_;          \
    int * const restrict * const restrict sim =                             \
            delete_sim ? allocate_int_table(2U,longestLen) : sim_;          \
    int * const restrict * const restrict len =                             \
            delete_len ? allocate_int_table(2U,longestLen) : len_;          \
    int * const restrict * const restrict del =                             \
            delete_del ? allocate_int_table(2U,longestLen) : del_;          \
    int * const restrict * const restrict ins =                             \
            delete_ins ? allocate_int_table(2U,longestLen) : ins_;          \
    cell_t maxCell = {NEG_INF, 0, 0, 0};   /* max cell in table */          \

#define AFFINE_GAP_DECL_SEMI                                                \
    cell_t lastCol = {NEG_INF, 0, 0, 0};   /* max cell in last col */       \
    cell_t lastRow = {NEG_INF, 0, 0, 0};   /* max cell in last row */       \

#define AFFINE_GAP_ASSERT       \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(scr);                \
    assert(mat);                \
    assert(sim);                \
    assert(len);                \
    assert(del);                \
    assert(ins);

#define AFFINE_GAP_INIT_CORNER                                              \
    scr[0][0] = 0;                                                          \
    mat[0][0] = 0;                                                          \
    sim[0][0] = 0;                                                          \
    len[0][0] = 0;                                                          \
    del[0][0] = NEG_INF;                                                    \
    ins[0][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY                              \
    for (j = 1; j <= s2Len; j++) {                                          \
        scr[0][j] = open + (j-1) * gap;                                     \
        mat[0][j] = 0;                                                      \
        sim[0][j] = 0;                                                      \
        len[0][j] = 0;                                                      \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY                           \
    for (j = 1; j <= s2Len; j++) {                                          \
        scr[0][j] = 0;                                                      \
        mat[0][j] = 0;                                                      \
        sim[0][j] = 0;                                                      \
        len[0][j] = 0;                                                      \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_BODY1                                                    \
    for (i = 1; i <= s1Len; ++i) {                                          \
        size_t cr = CROW(i);                                                \
        size_t pr = PROW(i);                                                \
        char ch1 = s1[i - 1];                                               \
        const int * const restrict BlosumRow = blosum_[MAP_BLOSUM_[ch1]];   \
                                                                            \
        int Nscore        = scr[pr][0];                                     \
        int Nmatches      = mat[pr][0];                                     \
        int Nsimilarities = sim[pr][0];                                     \
        int Nlength       = len[pr][0];                                     \
                                                                            \
        int Wscore        = open + (i-1) * gap;                             \
        int Wmatches      = 0;                                              \
        int Wsimilarities = 0;                                              \
        int Wlength       = 0;

#define AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY                              \
        /* init first column of 3 tables */                                 \
        scr[cr][0] = Wscore;                                                \
        mat[cr][0] = Wmatches;                                              \
        sim[cr][0] = Wsimilarities;                                         \
        len[cr][0] = Wlength;                                               \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY                           \
        /* init first column of 3 tables */                                 \
        scr[cr][0] = 0;                                                     \
        mat[cr][0] = Wmatches;                                              \
        sim[cr][0] = Wsimilarities;                                         \
        len[cr][0] = Wlength;                                               \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_BODY2(CALCULATE_SCORE)                                   \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int up = 0;                                                     \
            int left = 0;                                                   \
            int dig = 0;                                                    \
            char ch2 = s2[j - 1];                                           \
            int NWscore        = Nscore;                                    \
            int NWmatches      = Nmatches;                                  \
            int NWsimilarities = Nsimilarities;                             \
            int NWlength       = Nlength;                                   \
                                                                            \
            Nscore        = scr[pr][j];                                     \
            Nmatches      = mat[pr][j];                                     \
            Nsimilarities = sim[pr][j];                                     \
            Nlength       = len[pr][j];                                     \
                                                                            \
            tmp1 = CALCULATE_SCORE;                                         \
            up   = MAX(Nscore + open, del[pr][j]   + gap);                  \
            left = MAX(Wscore + open, ins[cr][j-1] + gap);                  \
            dig  = NWscore + tmp1;                                          \
                                                                            \
            del[cr][j] = up;                                                \
            ins[cr][j] = left;                                              \
                                                                            \
            if ((dig >= up) && (dig >= left)) {                             \
                Wscore = dig;                                               \
                Wmatches  = NWmatches + (ch1 == ch2);                       \
                Wsimilarities = NWsimilarities + (tmp1 > 0);                \
                Wlength  = NWlength + 1;                                    \
            } else if (up >= left) {                                        \
                Wscore = up;                                                \
                Wmatches  = Nmatches;                                       \
                Wsimilarities = Nsimilarities;                              \
                Wlength  = Nlength + 1;                                     \
            } else {                                                        \
                Wscore = left;                                              \
                Wmatches  = Wmatches;                                       \
                Wsimilarities  = Wsimilarities;                             \
                Wlength  = Wlength + 1;                                     \
            }                                                               \
                                                                            \
            scr[cr][j] = Wscore;                                            \
            mat[cr][j] = Wmatches;                                          \
            sim[cr][j] = Wsimilarities;                                     \
            len[cr][j] = Wlength;

#define AFFINE_GAP_BODY2_LOCAL                                              \
            if (scr[cr][j] > maxCell.score) {                               \
                maxCell.score        = scr[cr][j];                          \
                maxCell.matches      = mat[cr][j];                          \
                maxCell.similarities = sim[cr][j];                          \
                maxCell.length       = len[cr][j];                          \
            }                                                               \
            if (scr[cr][j] <= 0) {                                          \
                scr[cr][j] = Wscore = 0;                                    \
                mat[cr][j] = Wmatches = 0;                                  \
                sim[cr][j] = Wsimilarities = 0;                             \
                len[cr][j] = Wlength = 0;                                   \
            }

#define AFFINE_GAP_END_GLOBAL                                               \
        } /* end of j loop */                                               \
    } /* end of i loop */                                                   \
    maxCell.score        = scr[CROW(s1Len)][s2Len];                         \
    maxCell.matches      = mat[CROW(s1Len)][s2Len];                         \
    maxCell.similarities = sim[CROW(s1Len)][s2Len];                         \
    maxCell.length       = len[CROW(s1Len)][s2Len];

#define AFFINE_GAP_END_SEMI                                                 \
            /* track the maximum of last row */                             \
            if (i == s1Len && scr[cr][j] > lastRow.score) {                 \
                lastRow.score        = scr[cr][j];                          \
                lastRow.matches      = mat[cr][j];                          \
                lastRow.similarities = sim[cr][j];                          \
                lastRow.length       = len[cr][j];                          \
            }                                                               \
        } /* end of j loop */                                               \
        /* update the maximum of last column */                             \
        if (scr[cr][s2Len] > lastCol.score) {                               \
            lastCol.score        = scr[cr][s2Len];                          \
            lastCol.matches      = mat[cr][s2Len];                          \
            lastCol.similarities = sim[cr][s2Len];                          \
            lastCol.length       = len[cr][s2Len];                          \
        }                                                                   \
    } /* end of i loop */                                                   \
    maxCell = (lastCol.score > lastRow.score) ? lastCol : lastRow;

#define AFFINE_GAP_END_LOCAL                                                \
        } /* end of j loop */                                               \
    } /* end of i loop */

#define AFFINE_GAP_FREE                                                     \
    if (delete_scr) free_table(scr,2U);                                     \
    if (delete_mat) free_table(mat,2U);                                     \
    if (delete_sim) free_table(sim,2U);                                     \
    if (delete_len) free_table(len,2U);                                     \
    if (delete_del) free_table(del,2U);                                     \
    if (delete_ins) free_table(ins,2U);

#define AFFINE_GAP_RETURN                                                   \
    return maxCell;

#else /* USE_SIMILARITIES */

#define AFFINE_GAP_DECL                                                     \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    bool delete_scr = (scr_ == NULL);                                       \
    bool delete_mat = (mat_ == NULL);                                       \
    bool delete_len = (len_ == NULL);                                       \
    bool delete_del = (del_ == NULL);                                       \
    bool delete_ins = (ins_ == NULL);                                       \
    size_t longestLen = MAX(s1Len,s2Len);                                   \
    int * const restrict * const restrict scr =                             \
            delete_scr ? allocate_int_table(2U,longestLen) : scr_;          \
    int * const restrict * const restrict mat =                             \
            delete_mat ? allocate_int_table(2U,longestLen) : mat_;          \
    int * const restrict * const restrict len =                             \
            delete_len ? allocate_int_table(2U,longestLen) : len_;          \
    int * const restrict * const restrict del =                             \
            delete_del ? allocate_int_table(2U,longestLen) : del_;          \
    int * const restrict * const restrict ins =                             \
            delete_ins ? allocate_int_table(2U,longestLen) : ins_;          \
    cell_t maxCell = {NEG_INF, 0, 0};   /* max cell in table */             \

#define AFFINE_GAP_DECL_SEMI                                                \
    cell_t lastCol = {NEG_INF, 0, 0};   /* max cell in last col */          \
    cell_t lastRow = {NEG_INF, 0, 0};   /* max cell in last row */          \

#define AFFINE_GAP_ASSERT       \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(scr);                \
    assert(mat);                \
    assert(len);                \
    assert(del);                \
    assert(ins);

#define AFFINE_GAP_INIT_CORNER                                              \
    scr[0][0] = 0;                                                          \
    mat[0][0] = 0;                                                          \
    len[0][0] = 0;                                                          \
    del[0][0] = NEG_INF;                                                    \
    ins[0][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY                              \
    for (j = 1; j <= s2Len; j++) {                                          \
        scr[0][j] = open + (j-1) * gap;                                     \
        mat[0][j] = 0;                                                      \
        len[0][j] = 0;                                                      \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY                           \
    for (j = 1; j <= s2Len; j++) {                                          \
        scr[0][j] = 0;                                                      \
        mat[0][j] = 0;                                                      \
        len[0][j] = 0;                                                      \
        del[0][j] = NEG_INF;                                                \
        ins[0][j] = open + (j-1) * gap;                                     \
    }

#define AFFINE_GAP_BODY1                                                    \
    for (i = 1; i <= s1Len; ++i) {                                          \
        size_t cr = CROW(i);                                                \
        size_t pr = PROW(i);                                                \
        char ch1 = s1[i - 1];                                               \
        const int * const restrict BlosumRow = blosum_[MAP_BLOSUM_[ch1]];   \
                                                                            \
        int Nscore        = scr[pr][0];                                     \
        int Nmatches      = mat[pr][0];                                     \
        int Nlength       = len[pr][0];                                     \
                                                                            \
        int Wscore        = open + (i-1) * gap;                             \
        int Wmatches      = 0;                                              \
        int Wlength       = 0;

#define AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY                              \
        /* init first column of 3 tables */                                 \
        scr[cr][0] = Wscore;                                                \
        mat[cr][0] = Wmatches;                                              \
        len[cr][0] = Wlength;                                               \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY                           \
        /* init first column of 3 tables */                                 \
        scr[cr][0] = 0;                                                     \
        mat[cr][0] = Wmatches;                                              \
        len[cr][0] = Wlength;                                               \
        del[cr][0] = open + (i-1) * gap;                                    \
        ins[cr][0] = NEG_INF;

#define AFFINE_GAP_BODY2(CALCULATE_SCORE)                                   \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int up = 0;                                                     \
            int left = 0;                                                   \
            int dig = 0;                                                    \
            char ch2 = s2[j - 1];                                           \
            int NWscore        = Nscore;                                    \
            int NWmatches      = Nmatches;                                  \
            int NWlength       = Nlength;                                   \
                                                                            \
            Nscore        = scr[pr][j];                                     \
            Nmatches      = mat[pr][j];                                     \
            Nlength       = len[pr][j];                                     \
                                                                            \
            tmp1 = CALCULATE_SCORE;                                         \
            up   = MAX(Nscore + open, del[pr][j]   + gap);                  \
            left = MAX(Wscore + open, ins[cr][j-1] + gap);                  \
            dig  = NWscore + tmp1;                                          \
                                                                            \
            del[cr][j] = up;                                                \
            ins[cr][j] = left;                                              \
                                                                            \
            if ((dig >= up) && (dig >= left)) {                             \
                Wscore = dig;                                               \
                Wmatches  = NWmatches + (ch1 == ch2);                       \
                Wlength  = NWlength + 1;                                    \
            } else if (up >= left) {                                        \
                Wscore = up;                                                \
                Wmatches  = Nmatches;                                       \
                Wlength  = Nlength + 1;                                     \
            } else {                                                        \
                Wscore = left;                                              \
                Wmatches  = Wmatches;                                       \
                Wlength  = Wlength + 1;                                     \
            }                                                               \
                                                                            \
            scr[cr][j] = Wscore;                                            \
            mat[cr][j] = Wmatches;                                          \
            len[cr][j] = Wlength;

#define AFFINE_GAP_BODY2_LOCAL                                              \
            if (scr[cr][j] > maxCell.score) {                               \
                maxCell.score        = scr[cr][j];                          \
                maxCell.matches      = mat[cr][j];                          \
                maxCell.length       = len[cr][j];                          \
            }                                                               \
            if (scr[cr][j] <= 0) {                                          \
                scr[cr][j] = Wscore = 0;                                    \
                mat[cr][j] = Wmatches = 0;                                  \
                len[cr][j] = Wlength = 0;                                   \
            }

#define AFFINE_GAP_END_GLOBAL                                               \
        } /* end of j loop */                                               \
    } /* end of i loop */                                                   \
    maxCell.score        = scr[CROW(s1Len)][s2Len];                         \
    maxCell.matches      = mat[CROW(s1Len)][s2Len];                         \
    maxCell.length       = len[CROW(s1Len)][s2Len];

#define AFFINE_GAP_END_SEMI                                                 \
            /* track the maximum of last row */                             \
            if (i == s1Len && scr[cr][j] > lastRow.score) {                 \
                lastRow.score        = scr[cr][j];                          \
                lastRow.matches      = mat[cr][j];                          \
                lastRow.length       = len[cr][j];                          \
            }                                                               \
        } /* end of j loop */                                               \
        /* update the maximum of last column */                             \
        if (scr[cr][s2Len] > lastCol.score) {                               \
            lastCol.score        = scr[cr][s2Len];                          \
            lastCol.matches      = mat[cr][s2Len];                          \
            lastCol.length       = len[cr][s2Len];                          \
        }                                                                   \
    } /* end of i loop */                                                   \
    maxCell = (lastCol.score > lastRow.score) ? lastCol : lastRow;

#define AFFINE_GAP_END_LOCAL                                                \
        } /* end of j loop */                                               \
    } /* end of i loop */

#define AFFINE_GAP_FREE                                                     \
    if (delete_scr) free_table(scr,2U);                                     \
    if (delete_mat) free_table(mat,2U);                                     \
    if (delete_len) free_table(len,2U);                                     \
    if (delete_del) free_table(del,2U);                                     \
    if (delete_ins) free_table(ins,2U);

#define AFFINE_GAP_RETURN                                                   \
    return maxCell;

#endif /* USE_SIMILARITIES */



cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_global_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITH_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITH_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_END_GLOBAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_DECL_SEMI
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_END_SEMI
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        match_t callback,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(callback(ch1, ch2))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    assert(sub);
    assert(map);
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(SCORE(sub, map, first, ch1, ch2))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int match, int mismatch,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2((ch1 == ch2 ? match : mismatch))
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}

cell_t align_local_affine(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        int open, int gap,
        int * const restrict * const restrict scr_,
        int * const restrict * const restrict mat_,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim_,
#endif
        int * const restrict * const restrict len_,
        int * const restrict * const restrict del_,
        int * const restrict * const restrict ins_)
{
    AFFINE_GAP_DECL
    AFFINE_GAP_ASSERT
    AFFINE_GAP_INIT_CORNER
    AFFINE_GAP_INIT_FIRST_ROW_WITHOUT_PENALTY
    AFFINE_GAP_BODY1
    AFFINE_GAP_INIT_FIRST_COL_WITHOUT_PENALTY
    AFFINE_GAP_BODY2(BlosumRow[MAP_BLOSUM_[ch2]])
    AFFINE_GAP_BODY2_LOCAL
    AFFINE_GAP_END_LOCAL
    AFFINE_GAP_FREE
    AFFINE_GAP_RETURN
}


#define IS_EDGE_ASSERT  \
    assert(s1);         \
    assert(s2);         \
    assert(s1_len);     \
    assert(s2_len);     \

#define IS_EDGE_BODY(SELF_SCORE1, SELF_SCORE2)                              \
    if (s1_len > s2_len) {                                                  \
        max_len = s1_len;                                                   \
        self_score_ = SELF_SCORE1;                                          \
    }                                                                       \
    else {                                                                  \
        max_len = s2_len;                                                   \
        self_score_ = SELF_SCORE2;                                          \
    }                                                                       \
                                                                            \
    if (result.score <= 0) {                                                \
        return false;                                                       \
    }                                                                       \
    else if ((result.length * 100 >= AOL * int(max_len))                    \
             && (result.matches * 100 >= SIM * result.length)               \
             && (result.score * 100 >= OS * self_score_)) {                 \
        return true;                                                        \
    }                                                                       \
    else {                                                                  \
        return false;                                                       \
    }


bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score_,
        size_t &max_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first)
{
    IS_EDGE_ASSERT
    assert(sub);
    assert(map);

    IS_EDGE_BODY(self_score(s1, s1_len, sub, map, first),
                 self_score(s2, s2_len, sub, map, first))
}


bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score_,
        size_t &max_len,
        match_t callback)
{
    IS_EDGE_ASSERT

    IS_EDGE_BODY(self_score(s1, s1_len, callback),
                 self_score(s2, s2_len, callback))
}


bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score_,
        size_t &max_len,
        int match)
{
    IS_EDGE_ASSERT

    IS_EDGE_BODY(self_score(s1, s1_len, match),
                 self_score(s2, s2_len, match))
}


bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score_,
        size_t &max_len)
{
    IS_EDGE_ASSERT

    IS_EDGE_BODY(self_score_blosum(s1, s1_len),
                 self_score_blosum(s2, s2_len))
}


cell_t align_local_affine_ssw(
#if NOSWAP
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
#else
        const char * const restrict s2, size_t s2_len,
        const char * const restrict s1, size_t s1_len,
#endif
        int open, int gap)
{
    cell_t ret;
#if HAVE_EMMINTRIN_H
    s_profile *profile = NULL;
    int8_t *s1_num = new int8_t[s1_len];
    int8_t *s2_num = new int8_t[s2_len];
    s_align *result = NULL;

    /* This table is used to transform amino acid letters into numbers. */
    static const int8_t table[128] = {
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
    };

    /* initialize score matrix */
    for (size_t m = 0; m < s1_len; ++m) s1_num[m] = table[(int)s1[m]];
    for (size_t m = 0; m < s2_len; ++m) s2_num[m] = table[(int)s2[m]];
    profile = ssw_init(s1_num, s1_len, blosum__, 24, 2);
    result = ssw_align(profile, s2_num, s2_len, -open, -gap, 2, 0, 0, s1_len/2);

    ret.score = result->score1;
    ret.matches = 0;
#if USE_SIMILARITIES
    ret.similarities = 0;
#endif
    ret.length = 0;

    s_align* a = result;
    const char * const restrict ref_seq = s2;
    const char * const restrict read_seq = s1;

    if (a->cigar) {
        int32_t i, c = 0, qb = a->ref_begin1, pb = a->read_begin1;
        int32_t q = qb;
        int32_t p = pb;
        for (c = 0; c < a->cigarLen; ++c) {
            char letter = cigar_int_to_op(a->cigar[c]);
            uint32_t length = cigar_int_to_len(a->cigar[c]);
            for (i = 0; i < length; ++i){ 
                if (letter == 'M') {
                    int t1 = (int)*(ref_seq + q);
                    int t2 = (int)*(read_seq + p);
                    if (table[t1] == table[t2]) {
                        ret.matches += 1;
#if USE_SIMILARITIES
                        ret.similarities += 1;
#endif
                    }
                    else if (blosum__[table[t1]*24+table[t2]] > 0) {
#if USE_SIMILARITIES
                        ret.similarities += 1;
#endif
                    }
                    ++q;
                    ++p;
                } else {
                    if (letter == 'I') ++p;
                    else ++q;
                }
                ret.length += 1;
            }
        }
    }

    init_destroy(profile);
    if (result->cigar) free(result->cigar);
    free(result);
    delete [] s1_num;
    delete [] s2_num;
#else
    assert(0);
#endif

    return ret;
}


cell_t align_global_affine_fasta(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open, int gap)
{
    cell_t ret;
#if HAVE_EMMINTRIN_H
    int bias_signed = 0;
    unsigned short bias = 0;
    int32_t seg_len = 0;
    __m128i* profile_byte = NULL;
    int8_t* t = NULL;
    int32_t nt = 0;
    int32_t i = 0;
    int32_t j = 0;
    int32_t seg_num = 0;
    unsigned char *s1_num = new unsigned char[s1_len];
    unsigned char *s2_num = new unsigned char[s2_len];
    int score = 0;
    unsigned short ceiling = 255;
    struct f_struct f_str;

    /* This table is used to transform amino acid letters into numbers. */
    static const unsigned char table[128] = {
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
    };

    /* initialize score matrix */
    for (size_t m = 0; m < s1_len; ++m) s1_num[m] = table[(int)s1[m]];
    for (size_t m = 0; m < s2_len; ++m) s2_num[m] = table[(int)s2[m]];
    /* find the bias; limit is 24x24 for blosum matrix */
    for (i=0; i<24*24; ++i) {
        if (blosum__[i] < bias_signed) {
            bias_signed = blosum__[i];
        }
    }
    bias = abs(bias_signed);

    /* create query profile for byte routine */
    {
        /* Split the 128 bit register into 16 pieces. Each piece is 8
         * bit. Split the read into 16 segments. Calculat 16 segments in
         * parallel.  */
        seg_len = (s1_len + 15) / 16;
        profile_byte = new __m128i[24 * seg_len];
        t = (int8_t*)profile_byte; 

        /* Generate query profile rearrange query sequence & calculate
         * the weight of match/mismatch */
        for (nt = 0; nt < 24; ++nt) {
            for (i = 0; i < seg_len; ++i) {
                j = i;
                for (seg_num = 0; seg_num < 16; ++seg_num) {
                    *t++ = j>= s1_len ? bias : blosum__[nt * 24 + s1_num[j]] + bias;
                    j += seg_len;
                }
            }
        }
    }

    /* create workspace for table */
    {
        seg_len = (s1_len + 7) / 8;
        f_str.workspace = new __m128i[3*seg_len];
    }

    score = global_sse2_byte(
            s1_len, (unsigned char*)profile_byte, s2_num, s2_len,
            -open, -gap, ceiling, bias, &f_str);
    ret.score = score;
#else
    assert(0);
#endif
    return ret;
}

#define BLOSUM0_(ch1, ch2) (matrow0[(ch2)])
#define BLOSUM1_(ch1, ch2) (matrow1[(ch2)])
#define BLOSUM2_(ch1, ch2) (matrow2[(ch2)])
#define BLOSUM3_(ch1, ch2) (matrow3[(ch2)])
#define BLOSUM4_(ch1, ch2) (matrow4[(ch2)])
#define BLOSUM5_(ch1, ch2) (matrow5[(ch2)])
#define BLOSUM6_(ch1, ch2) (matrow6[(ch2)])
#define BLOSUM7_(ch1, ch2) (matrow7[(ch2)])

/* on OSX, _mm_extract_epi16 got wrong answer, but taking the union of
 * the vector and extracting that way seemed to work... */
#define EXTRACT extract
//#define EXTRACT _mm_extract_epi16

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
                    (m) = EXTRACT((vm), 0)

static inline int extract(const __m128i &m, const int &pos)
{
    union {
        __m128i m;
        int16_t v[8];
    } tmp;
    tmp.m = m;
    return tmp.v[pos];
}


/* shift given vector v and insert val */
static inline __m128i vshift16(const __m128i &v)
{
    return _mm_srli_si128(v, 2);
}


/* shift given vector v and insert val */
static inline __m128i vshift16(const __m128i &v, int val)
{
    __m128i ret = _mm_srli_si128(v, 2);
    ret = _mm_insert_epi16(ret, val, 7);
    return ret;
}


/* shift given vector v and insert val */
static inline int vextract16(const __m128i &v, int offset)
{
    switch (offset) {
        case 0: return EXTRACT(v, 0);
        case 1: return EXTRACT(v, 1);
        case 2: return EXTRACT(v, 2);
        case 3: return EXTRACT(v, 3);
        case 4: return EXTRACT(v, 4);
        case 5: return EXTRACT(v, 5);
        case 6: return EXTRACT(v, 6);
        case 7: return EXTRACT(v, 7);
        default: assert(0);
    }
}


cell_t align_global_affine_sse(
        const char * const restrict _s1, const size_t _s1Len,
        const char * const restrict _s2, const size_t _s2Len,
        const int _open, const int _gap,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    cell_t result;
    int s1Len = int(_s1Len);
    int s2Len = int(_s2Len);
    const int open = -_open;
    const int gap = -_gap;
    int score = 0;
    int match = 0;
    int length = 0;
    int * const restrict s1 = new int[s1Len+8];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int i=s1Len; i<s1Len+8; ++i) {
        s1[i] = 23;
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
    mch_pr[7] = 0;
    len_pr[7] = 0;
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = -open -(j-8)*gap;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (int i=1; i<=s1Len; i+=8) {
        int j;
        __m128i NWscore  = _mm_set1_epi16(NEG_INF);
        __m128i NWmatch  = _mm_set1_epi16(NEG_INF);
        __m128i NWlength = _mm_set1_epi16(NEG_INF);
        __m128i Nscore   = _mm_set1_epi16(NEG_INF);
        __m128i Nmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Nlength  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore   = _mm_set1_epi16(NEG_INF);
        __m128i Wmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Wlength  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl     = _mm_set1_epi16(NEG_INF);
        __m128i vDel     = _mm_set1_epi16(NEG_INF);
        __m128i vIns     = _mm_set1_epi16(NEG_INF);
        __m128i vMat;
        __m128i vs1;
        __m128i vs2;
        __m128i vOne     = _mm_set1_epi16(1);
        __m128i case1not;
        __m128i case2not;
        __m128i case2;
        __m128i case3;
        __m128i Cscore;
        __m128i Cmatch;
        __m128i Clength;

        const int * const restrict matrow0 = blosum_[s1[i-1+0]];
        const int * const restrict matrow1 = blosum_[s1[i-1+1]];
        const int * const restrict matrow2 = blosum_[s1[i-1+2]];
        const int * const restrict matrow3 = blosum_[s1[i-1+3]];
        const int * const restrict matrow4 = blosum_[s1[i-1+4]];
        const int * const restrict matrow5 = blosum_[s1[i-1+5]];
        const int * const restrict matrow6 = blosum_[s1[i-1+6]];
        const int * const restrict matrow7 = blosum_[s1[i-1+7]];

        vs1 = _mm_set_epi16(
                s1[i-1+0],
                s1[i-1+1],
                s1[i-1+2],
                s1[i-1+3],
                s1[i-1+4],
                s1[i-1+5],
                s1[i-1+6],
                s1[i-1+7]
                );

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Nmatch = _mm_insert_epi16(Nscore, mch_pr[j+7], 7);
        Nlength= _mm_insert_epi16(Nscore, len_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, -open-(i-1)*gap, 7);
        Wmatch = _mm_insert_epi16(Wscore, 0, 7);
        Wlength= _mm_insert_epi16(Wscore, 0, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);

        /* j = 1 */
        j = 1;
        /* this block never changes, so let's not copy-and-paste */
#ifdef SETUP_BLOCK
#undef SETUP_BLOCK
#endif
#define SETUP_BLOCK                                 \
        NWscore = Nscore;                           \
        NWmatch = Nmatch;                           \
        NWlength= Nlength;                          \
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);    \
        Nmatch  = vshift16(Wmatch, mch_pr[j+7]);    \
        Nlength = vshift16(Wlength,len_pr[j+7]);    \
        vDel    = vshift16(vDel,   del_pr[j+7]);    \
        vDel = _mm_max_epi16(                       \
                _mm_sub_epi16(Nscore, vOpen),       \
                _mm_sub_epi16(vDel, vGap));         \
        vIns = _mm_max_epi16(                       \
                _mm_sub_epi16(Wscore,vOpen),        \
                _mm_sub_epi16(vIns,vGap));
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1],s2[j-1]),
            0,
            0,
            0,
            0,
            0,
            0,
            0
        );
        /* we reuse this logic throughout this function,
         * best to not copy-and-paste everywhere */
#ifdef CONDITIONAL_BLOCK
#undef CONDITIONAL_BLOCK
#endif
#define CONDITIONAL_BLOCK                                                   \
        vTbl = _mm_add_epi16(NWscore, vMat);                                \
        case1not = _mm_or_si128(                                            \
                _mm_cmplt_epi16(vTbl,vDel),_mm_cmplt_epi16(vTbl,vIns));     \
        case2not = _mm_cmplt_epi16(vDel,vIns);                              \
        case2 = _mm_andnot_si128(case2not,case1not);                        \
        case3 = _mm_and_si128(case1not,case2not);                           \
        Cscore = _mm_andnot_si128(case1not, vTbl);                          \
        Cmatch = _mm_andnot_si128(case1not,                                 \
                    _mm_add_epi16(NWmatch, _mm_and_si128(                   \
                            _mm_cmpeq_epi16(vs1,vs2),vOne)));               \
        Clength= _mm_andnot_si128(case1not, _mm_add_epi16(NWlength, vOne)); \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case2, vDel));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case2,                  \
                    _mm_add_epi16(Nlength, vOne)));                         \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case3, vIns));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case3,                  \
                    _mm_add_epi16(Wlength, vOne)));                         \
        Wscore = vTbl = Cscore;                                             \
        Wmatch = Cmatch;                                                    \
        Wlength= Clength;
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+0)*gap, 6);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 6);
        Wlength= _mm_insert_epi16(Wlength,0, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);

        /* j = 2 */
        j = 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            0,
            0,
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+1)*gap, 5);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 5);
        Wlength= _mm_insert_epi16(Wlength,0, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);

        /* j = 3 */
        j = 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            0,
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+2)*gap, 4);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 4);
        Wlength= _mm_insert_epi16(Wlength,0, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);

        /* j = 4 */
        j = 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+3)*gap, 3);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 3);
        Wlength= _mm_insert_epi16(Wlength,0, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);

        /* j = 5 */
        j = 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+4)*gap, 2);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 2);
        Wlength= _mm_insert_epi16(Wlength,0, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);

        /* j = 6 */
        j = 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+5)*gap, 1);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 1);
        Wlength= _mm_insert_epi16(Wlength,0, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);

        /* j = 7 */
        j = 7;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            BLOSUM6_(s1[i-1+6],s2[j-1-6]),
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+6)*gap, 0);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 0);
        Wlength= _mm_insert_epi16(Wlength,0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = -open-(i+6)*gap;
        mch_pr[7] = 0;
        len_pr[7] = 0;

        for (j=8; j<=s2Len; ++j) {
            SETUP_BLOCK
            vs2 = vshift16(vs2, s2[j-1]);
            vMat = _mm_set_epi16(
                    BLOSUM0_(s1[i-1+0],s2[j-1-0]),
                    BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                    BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                    BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                    BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                    BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                    BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                    BLOSUM7_(s1[i-1+7],s2[j-1-7])
                    );
            CONDITIONAL_BLOCK
            tbl_pr[j] = EXTRACT(vTbl, 0);
            mch_pr[j] = EXTRACT(Wmatch, 0);
            len_pr[j] = EXTRACT(Wlength, 0);
            del_pr[j] = EXTRACT(vDel, 0);
        }
        if (i+0 == s1Len) {
            score = EXTRACT(vTbl, 7);
            match = EXTRACT(Wmatch, 7);
            length= EXTRACT(Wlength, 7);
        }

        /* j = s2Len + 1 */
        j = s2Len + 1;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+1 == s1Len) {
            score = EXTRACT(vTbl, 6);
            match = EXTRACT(Wmatch, 6);
            length= EXTRACT(Wlength, 6);
        }

        /* j = s2Len + 2 */
        j = s2Len + 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+2 == s1Len) {
            score = EXTRACT(vTbl, 5);
            match = EXTRACT(Wmatch, 5);
            length= EXTRACT(Wlength, 5);
        }

        /* j = s2Len + 3 */
        j = s2Len + 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+3 == s1Len) {
            score = EXTRACT(vTbl, 4);
            match = EXTRACT(Wmatch, 4);
            length= EXTRACT(Wlength, 4);
        }

        /* j = s2Len + 4 */
        j = s2Len + 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+4 == s1Len) {
            score = EXTRACT(vTbl, 3);
            match = EXTRACT(Wmatch, 3);
            length= EXTRACT(Wlength, 3);
        }

        /* j = s2Len + 5 */
        j = s2Len + 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+5 == s1Len) {
            score = EXTRACT(vTbl, 2);
            match = EXTRACT(Wmatch, 2);
            length= EXTRACT(Wlength, 2);
        }

        /* j = s2Len + 6 */
        j = s2Len + 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+6 == s1Len) {
            score = EXTRACT(vTbl, 1);
            match = EXTRACT(Wmatch, 1);
            length= EXTRACT(Wlength, 1);
        }

        /* j = s2Len + 7 */
        j = s2Len + 7;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (i+7 == s1Len) {
            score = EXTRACT(vTbl, 0);
            match = EXTRACT(Wmatch, 0);
            length= EXTRACT(Wlength, 0);
        }
    }

    delete [] s1;
    delete [] s2;
    result.score = score;
    result.matches = match;
    result.length = length;
#if USE_SIMILARITIES
    result.similarities = 0;
#endif
    return result;
}


cell_t align_semi_affine_sse(
        const char * const restrict _s1, const size_t _s1Len,
        const char * const restrict _s2, const size_t _s2Len,
        const int _open, const int _gap,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    cell_t result;
    int s1Len = int(_s1Len);
    int s2Len = int(_s2Len);
    const int open = -_open;
    const int gap = -_gap;
    int score = NEG_INF;
    int match = 0;
    int length = 0;
    int * const restrict s1 = new int[s1Len+8];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int i=s1Len; i<s1Len+8; ++i) {
        s1[i] = 23;
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
    mch_pr[7] = 0;
    len_pr[7] = 0;
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (int i=1; i<=s1Len; i+=8) {
        bool last_pass = (i+8>s1Len);
        int offset = 7 - (s1Len - i);
        int j;
        __m128i NWscore  = _mm_set1_epi16(NEG_INF);
        __m128i NWmatch  = _mm_set1_epi16(NEG_INF);
        __m128i NWlength = _mm_set1_epi16(NEG_INF);
        __m128i Nscore   = _mm_set1_epi16(NEG_INF);
        __m128i Nmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Nlength  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore   = _mm_set1_epi16(NEG_INF);
        __m128i Wmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Wlength  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl     = _mm_set1_epi16(NEG_INF);
        __m128i vDel     = _mm_set1_epi16(NEG_INF);
        __m128i vIns     = _mm_set1_epi16(NEG_INF);
        __m128i vMat;
        __m128i vs1;
        __m128i vs2;
        __m128i vOne     = _mm_set1_epi16(1);
        __m128i case1not;
        __m128i case2not;
        __m128i case2;
        __m128i case3;
        __m128i Cscore;
        __m128i Cmatch;
        __m128i Clength;

        const int * const restrict matrow0 = blosum_[s1[i-1+0]];
        const int * const restrict matrow1 = blosum_[s1[i-1+1]];
        const int * const restrict matrow2 = blosum_[s1[i-1+2]];
        const int * const restrict matrow3 = blosum_[s1[i-1+3]];
        const int * const restrict matrow4 = blosum_[s1[i-1+4]];
        const int * const restrict matrow5 = blosum_[s1[i-1+5]];
        const int * const restrict matrow6 = blosum_[s1[i-1+6]];
        const int * const restrict matrow7 = blosum_[s1[i-1+7]];

        vs1 = _mm_set_epi16(
                s1[i-1+0],
                s1[i-1+1],
                s1[i-1+2],
                s1[i-1+3],
                s1[i-1+4],
                s1[i-1+5],
                s1[i-1+6],
                s1[i-1+7]
                );

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Nmatch = _mm_insert_epi16(Nscore, mch_pr[j+7], 7);
        Nlength= _mm_insert_epi16(Nscore, len_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, 0, 7);
        Wmatch = _mm_insert_epi16(Wscore, 0, 7);
        Wlength= _mm_insert_epi16(Wscore, 0, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);

        /* j = 1 */
        j = 1;
        /* this block never changes, so let's not copy-and-paste */
#ifdef SETUP_BLOCK
#undef SETUP_BLOCK
#endif
#define SETUP_BLOCK                                 \
        NWscore = Nscore;                           \
        NWmatch = Nmatch;                           \
        NWlength= Nlength;                          \
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);    \
        Nmatch  = vshift16(Wmatch, mch_pr[j+7]);    \
        Nlength = vshift16(Wlength,len_pr[j+7]);    \
        vDel    = vshift16(vDel,   del_pr[j+7]);    \
        vDel = _mm_max_epi16(                       \
                _mm_sub_epi16(Nscore, vOpen),       \
                _mm_sub_epi16(vDel, vGap));         \
        vIns = _mm_max_epi16(                       \
                _mm_sub_epi16(Wscore,vOpen),        \
                _mm_sub_epi16(vIns,vGap));
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1],s2[j-1]),
            0,
            0,
            0,
            0,
            0,
            0,
            0
        );
        /* we reuse this logic throughout this function,
         * best to not copy-and-paste everywhere */
#ifdef CONDITIONAL_BLOCK
#undef CONDITIONAL_BLOCK
#endif
#define CONDITIONAL_BLOCK                                                   \
        vTbl = _mm_add_epi16(NWscore, vMat);                                \
        case1not = _mm_or_si128(                                            \
                _mm_cmplt_epi16(vTbl,vDel),_mm_cmplt_epi16(vTbl,vIns));     \
        case2not = _mm_cmplt_epi16(vDel,vIns);                              \
        case2 = _mm_andnot_si128(case2not,case1not);                        \
        case3 = _mm_and_si128(case1not,case2not);                           \
        Cscore = _mm_andnot_si128(case1not, vTbl);                          \
        Cmatch = _mm_andnot_si128(case1not,                                 \
                    _mm_add_epi16(NWmatch, _mm_and_si128(                   \
                            _mm_cmpeq_epi16(vs1,vs2),vOne)));               \
        Clength= _mm_andnot_si128(case1not, _mm_add_epi16(NWlength, vOne)); \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case2, vDel));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case2,                  \
                    _mm_add_epi16(Nlength, vOne)));                         \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case3, vIns));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case3,                  \
                    _mm_add_epi16(Wlength, vOne)));                         \
        Wscore = vTbl = Cscore;                                             \
        Wmatch = Cmatch;                                                    \
        Wlength= Clength;
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 6);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 6);
        Wlength= _mm_insert_epi16(Wlength,0, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        /* j = 2 */
        j = 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            0,
            0,
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 5);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 5);
        Wlength= _mm_insert_epi16(Wlength,0, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        /* j = 3 */
        j = 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            0,
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 4);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 4);
        Wlength= _mm_insert_epi16(Wlength,0, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        /* j = 4 */
        j = 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 3);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 3);
        Wlength= _mm_insert_epi16(Wlength,0, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        /* j = 5 */
        j = 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 2);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 2);
        Wlength= _mm_insert_epi16(Wlength,0, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        /* j = 6 */
        j = 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 1);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 1);
        Wlength= _mm_insert_epi16(Wlength,0, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        /* j = 7 */
        j = 7;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            BLOSUM6_(s1[i-1+6],s2[j-1-6]),
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 0);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 0);
        Wlength= _mm_insert_epi16(Wlength,0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = 0;
        mch_pr[7] = 0;
        len_pr[7] = 0;
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 7);
                length= EXTRACT(Wlength, 7);
            }
        }

        for (j=8; j<=s2Len; ++j) {
            SETUP_BLOCK
            vs2 = vshift16(vs2, s2[j-1]);
            vMat = _mm_set_epi16(
                    BLOSUM0_(s1[i-1+0],s2[j-1-0]),
                    BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                    BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                    BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                    BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                    BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                    BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                    BLOSUM7_(s1[i-1+7],s2[j-1-7])
                    );
            CONDITIONAL_BLOCK
            tbl_pr[j] = EXTRACT(vTbl, 0);
            mch_pr[j] = EXTRACT(Wmatch, 0);
            len_pr[j] = EXTRACT(Wlength, 0);
            del_pr[j] = EXTRACT(vDel, 0);
            if (last_pass) {
                int tmp = vextract16(vTbl, offset);
                if (tmp > score) {
                    score = tmp;
                    match = vextract16(Wmatch, offset);
                    length= vextract16(Wlength, offset);
                }
            }
            if (j == s2Len) {
                int tmp = EXTRACT(vTbl, 7);
                if (tmp > score) {
                    score = tmp;
                    match = EXTRACT(Wmatch, 7);
                    length= EXTRACT(Wlength, 7);
                }
            }
        }

        /* j = s2Len + 1 */
        j = s2Len + 1;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 6);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 6);
                length= EXTRACT(Wlength, 6);
            }
        }

        /* j = s2Len + 2 */
        j = s2Len + 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 5);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 5);
                length= EXTRACT(Wlength, 5);
            }
        }

        /* j = s2Len + 3 */
        j = s2Len + 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 4);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 4);
                length= EXTRACT(Wlength, 4);
            }
        }

        /* j = s2Len + 4 */
        j = s2Len + 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 3);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 3);
                length= EXTRACT(Wlength, 3);
            }
        }

        /* j = s2Len + 5 */
        j = s2Len + 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 2);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 2);
                length= EXTRACT(Wlength, 2);
            }
        }

        /* j = s2Len + 6 */
        j = s2Len + 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 1);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 1);
                length= EXTRACT(Wlength, 1);
            }
        }

        /* j = s2Len + 7 */
        j = s2Len + 7;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            if (tmp > score) {
                score = tmp;
                match = vextract16(Wmatch, offset);
                length= vextract16(Wlength, offset);
            }
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 0);
            if (tmp > score) {
                score = tmp;
                match = EXTRACT(Wmatch, 0);
                length= EXTRACT(Wlength, 0);
            }
        }
    }

    delete [] s1;
    delete [] s2;
    result.score = score;
    result.matches = match;
    result.length = length;
#if USE_SIMILARITIES
    result.similarities = 0;
#endif
    return result;
}


cell_t align_local_affine_sse(
        const char * const restrict _s1, const size_t _s1Len,
        const char * const restrict _s2, const size_t _s2Len,
        const int _open, const int _gap,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    cell_t result;
    int s1Len = int(_s1Len);
    int s2Len = int(_s2Len);
    const int open = -_open;
    const int gap = -_gap;
    int score = NEG_INF;
    int match = 0;
    int length = 0;
    __m128i vScore = _mm_setzero_si128();
    __m128i vMatch = _mm_setzero_si128();
    __m128i vLength= _mm_setzero_si128();
    int * const restrict s1 = new int[s1Len+8];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int i=s1Len; i<s1Len+8; ++i) {
        s1[i] = 23;
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
    mch_pr[7] = 0;
    len_pr[7] = 0;
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (int i=1; i<=s1Len; i+=8) {
        bool last_pass = (i+8>=s1Len);
        int offset = 7 - (s1Len - i);
        int j;
        __m128i NWscore  = _mm_set1_epi16(NEG_INF);
        __m128i NWmatch  = _mm_set1_epi16(NEG_INF);
        __m128i NWlength = _mm_set1_epi16(NEG_INF);
        __m128i Nscore   = _mm_set1_epi16(NEG_INF);
        __m128i Nmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Nlength  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore   = _mm_set1_epi16(NEG_INF);
        __m128i Wmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Wlength  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl     = _mm_set1_epi16(NEG_INF);
        __m128i vDel     = _mm_set1_epi16(NEG_INF);
        __m128i vIns     = _mm_set1_epi16(NEG_INF);
        __m128i vMat;
        __m128i vs1;
        __m128i vs2;
        __m128i vOne     = _mm_set1_epi16(1);
        __m128i case1not;
        __m128i case2not;
        __m128i case2;
        __m128i case3;
        __m128i Cscore;
        __m128i Cmatch;
        __m128i Clength;
        __m128i mask;

        const int * const restrict matrow0 = blosum_[s1[i-1+0]];
        const int * const restrict matrow1 = blosum_[s1[i-1+1]];
        const int * const restrict matrow2 = blosum_[s1[i-1+2]];
        const int * const restrict matrow3 = blosum_[s1[i-1+3]];
        const int * const restrict matrow4 = blosum_[s1[i-1+4]];
        const int * const restrict matrow5 = blosum_[s1[i-1+5]];
        const int * const restrict matrow6 = blosum_[s1[i-1+6]];
        const int * const restrict matrow7 = blosum_[s1[i-1+7]];

        vs1 = _mm_set_epi16(
                s1[i-1+0],
                s1[i-1+1],
                s1[i-1+2],
                s1[i-1+3],
                s1[i-1+4],
                s1[i-1+5],
                s1[i-1+6],
                s1[i-1+7]
                );

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Nmatch = _mm_insert_epi16(Nscore, mch_pr[j+7], 7);
        Nlength= _mm_insert_epi16(Nscore, len_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, 0, 7);
        Wmatch = _mm_insert_epi16(Wscore, 0, 7);
        Wlength= _mm_insert_epi16(Wscore, 0, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);

        /* j = 1 */
        j = 1;
        /* this block never changes, so let's not copy-and-paste */
#ifdef SETUP_BLOCK
#undef SETUP_BLOCK
#endif
#define SETUP_BLOCK                                 \
        NWscore = Nscore;                           \
        NWmatch = Nmatch;                           \
        NWlength= Nlength;                          \
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);    \
        Nmatch  = vshift16(Wmatch, mch_pr[j+7]);    \
        Nlength = vshift16(Wlength,len_pr[j+7]);    \
        vDel    = vshift16(vDel,   del_pr[j+7]);    \
        vDel = _mm_max_epi16(                       \
                _mm_sub_epi16(Nscore, vOpen),       \
                _mm_sub_epi16(vDel, vGap));         \
        vIns = _mm_max_epi16(                       \
                _mm_sub_epi16(Wscore,vOpen),        \
                _mm_sub_epi16(vIns,vGap));
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1],s2[j-1]),
            0,
            0,
            0,
            0,
            0,
            0,
            0
        );
        /* we reuse this logic throughout this function,
         * best to not copy-and-paste everywhere */
#ifdef CONDITIONAL_BLOCK
#undef CONDITIONAL_BLOCK
#endif
#define CONDITIONAL_BLOCK                                                   \
        vTbl = _mm_add_epi16(NWscore, vMat);                                \
        case1not = _mm_or_si128(                                            \
                _mm_cmplt_epi16(vTbl,vDel),_mm_cmplt_epi16(vTbl,vIns));     \
        case2not = _mm_cmplt_epi16(vDel,vIns);                              \
        case2 = _mm_andnot_si128(case2not,case1not);                        \
        case3 = _mm_and_si128(case1not,case2not);                           \
        Cscore = _mm_andnot_si128(case1not, vTbl);                          \
        Cmatch = _mm_andnot_si128(case1not,                                 \
                    _mm_add_epi16(NWmatch, _mm_and_si128(                   \
                            _mm_cmpeq_epi16(vs1,vs2),vOne)));               \
        Clength= _mm_andnot_si128(case1not, _mm_add_epi16(NWlength, vOne)); \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case2, vDel));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case2,                  \
                    _mm_add_epi16(Nlength, vOne)));                         \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case3, vIns));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case3,                  \
                    _mm_add_epi16(Wlength, vOne)));                         \
        mask = _mm_cmplt_epi16(Cscore,vOne);                                \
        Cscore = _mm_andnot_si128(mask, Cscore);                            \
        Cmatch = _mm_andnot_si128(mask, Cmatch);                            \
        Clength= _mm_andnot_si128(mask, Clength);                           \
        mask = _mm_cmpgt_epi16(Cscore,vScore);                              \
        vScore = _mm_andnot_si128(mask, vScore);                            \
        vScore = _mm_or_si128(vScore, _mm_and_si128(mask, Cscore));         \
        vMatch = _mm_andnot_si128(mask, vMatch);                            \
        vMatch = _mm_or_si128(vMatch, _mm_and_si128(mask, Cmatch));         \
        vLength= _mm_andnot_si128(mask, vLength);                           \
        vLength= _mm_or_si128(vLength,_mm_and_si128(mask, Clength));        \
        Wscore = vTbl = Cscore;                                             \
        Wmatch = Cmatch;                                                    \
        Wlength= Clength;
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 6);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 6);
        Wlength= _mm_insert_epi16(Wlength,0, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);

        /* j = 2 */
        j = 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            0,
            0,
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 5);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 5);
        Wlength= _mm_insert_epi16(Wlength,0, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);

        /* j = 3 */
        j = 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            0,
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 4);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 4);
        Wlength= _mm_insert_epi16(Wlength,0, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);

        /* j = 4 */
        j = 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            0,
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 3);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 3);
        Wlength= _mm_insert_epi16(Wlength,0, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);

        /* j = 5 */
        j = 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            0,
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 2);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 2);
        Wlength= _mm_insert_epi16(Wlength,0, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);

        /* j = 6 */
        j = 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            0,
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 1);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 1);
        Wlength= _mm_insert_epi16(Wlength,0, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);

        /* j = 7 */
        j = 7;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            BLOSUM6_(s1[i-1+6],s2[j-1-6]),
            0
        );
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, 0, 0);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 0);
        Wlength= _mm_insert_epi16(Wlength,0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = 0;
        mch_pr[7] = 0;
        len_pr[7] = 0;

        for (j=8; j<=s2Len; ++j) {
            SETUP_BLOCK
            vs2 = vshift16(vs2, s2[j-1]);
            vMat = _mm_set_epi16(
                    BLOSUM0_(s1[i-1+0],s2[j-1-0]),
                    BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                    BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                    BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                    BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                    BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                    BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                    BLOSUM7_(s1[i-1+7],s2[j-1-7])
                    );
            CONDITIONAL_BLOCK
            tbl_pr[j] = EXTRACT(vTbl, 0);
            mch_pr[j] = EXTRACT(Wmatch, 0);
            len_pr[j] = EXTRACT(Wlength, 0);
            del_pr[j] = EXTRACT(vDel, 0);
        }

        /* j = s2Len + 1 */
        __m128i jMask = _mm_cmpeq_epi16(vOne,vOne);
        j = s2Len + 1;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        /* add jMask to mask out invalid j indices gt s2Len */
#ifdef CONDITIONAL_BLOCK
#undef CONDITIONAL_BLOCK
#endif
#define CONDITIONAL_BLOCK                                                   \
        vTbl = _mm_add_epi16(NWscore, vMat);                                \
        case1not = _mm_or_si128(                                            \
                _mm_cmplt_epi16(vTbl,vDel),_mm_cmplt_epi16(vTbl,vIns));     \
        case2not = _mm_cmplt_epi16(vDel,vIns);                              \
        case2 = _mm_andnot_si128(case2not,case1not);                        \
        case3 = _mm_and_si128(case1not,case2not);                           \
        Cscore = _mm_andnot_si128(case1not, vTbl);                          \
        Cmatch = _mm_andnot_si128(case1not,                                 \
                    _mm_add_epi16(NWmatch, _mm_and_si128(                   \
                            _mm_cmpeq_epi16(vs1,vs2),vOne)));               \
        Clength= _mm_andnot_si128(case1not, _mm_add_epi16(NWlength, vOne)); \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case2, vDel));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case2,                  \
                    _mm_add_epi16(Nlength, vOne)));                         \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case3, vIns));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case3,                  \
                    _mm_add_epi16(Wlength, vOne)));                         \
        jMask = _mm_srli_si128(jMask, 2);                                   \
        Cscore = _mm_and_si128(jMask, Cscore);                              \
        mask = _mm_cmplt_epi16(Cscore,vOne);                                \
        Cscore = _mm_andnot_si128(mask, Cscore);                            \
        Cmatch = _mm_andnot_si128(mask, Cmatch);                            \
        Clength= _mm_andnot_si128(mask, Clength);                           \
        mask = _mm_cmpgt_epi16(Cscore,vScore);                              \
        vScore = _mm_andnot_si128(mask, vScore);                            \
        vScore = _mm_or_si128(vScore, _mm_and_si128(mask, Cscore));         \
        vMatch = _mm_andnot_si128(mask, vMatch);                            \
        vMatch = _mm_or_si128(vMatch, _mm_and_si128(mask, Cmatch));         \
        vLength= _mm_andnot_si128(mask, vLength);                           \
        vLength= _mm_or_si128(vLength,_mm_and_si128(mask, Clength));        \
        Wscore = vTbl = Cscore;                                             \
        Wmatch = Cmatch;                                                    \
        Wlength= Clength;
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);

        /* j = s2Len + 2 */
        j = s2Len + 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);

        /* j = s2Len + 3 */
        j = s2Len + 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);

        /* j = s2Len + 4 */
        j = s2Len + 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);

        /* j = s2Len + 5 */
        j = s2Len + 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);

        /* j = s2Len + 6 */
        j = s2Len + 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);

        /* j = s2Len + 7 */
        j = s2Len + 7;
        SETUP_BLOCK
        vs2 = vshift16(vs2);
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
    }

    {
        union {
            __m128i m;
            int16_t v[8];
        } tScore, tMatch, tLength;
        tScore.m = vScore;
        tMatch.m = vMatch;
        tLength.m= vLength;
        score = tScore.v[7];
        match = tMatch.v[7];
        length= tLength.v[7];
        if (tScore.v[1] > score) {
            score = tScore.v[6];
            match = tMatch.v[6];
            length= tLength.v[6];
        }
        if (tScore.v[2] > score) {
            score = tScore.v[5];
            match = tMatch.v[5];
            length= tLength.v[5];
        }
        if (tScore.v[3] > score) {
            score = tScore.v[4];
            match = tMatch.v[4];
            length= tLength.v[4];
        }
        if (tScore.v[4] > score) {
            score = tScore.v[3];
            match = tMatch.v[3];
            length= tLength.v[3];
        }
        if (tScore.v[5] > score) {
            score = tScore.v[2];
            match = tMatch.v[2];
            length= tLength.v[2];
        }
        if (tScore.v[6] > score) {
            score = tScore.v[1];
            match = tMatch.v[1];
            length= tLength.v[1];
        }
        if (tScore.v[7] > score) {
            score = tScore.v[0];
            match = tMatch.v[0];
            length= tLength.v[0];
        }
    }
    delete [] s1;
    delete [] s2;
    result.score = score;
    result.matches = match;
    result.length = length;
#if USE_SIMILARITIES
    result.similarities = 0;
#endif
    return result;
}


}; /* namespace pgraph */

