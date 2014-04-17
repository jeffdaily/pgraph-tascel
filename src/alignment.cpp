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
static const int NEG_INF = INT_MIN / 2;


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
    ret.length = 0;

    s_align* a = result;
    const char * const restrict ref_seq = s2;
    const char * const restrict read_seq = s1;

    if (a->cigar) {
        int32_t i, c = 0, qb = a->ref_begin1, pb = a->read_begin1;
        int32_t q = qb;
        int32_t p = pb;
        for (c = 0; c < a->cigarLen; ++c) {
            int32_t letter = 0xf&*(a->cigar + c);
            int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
            for (i = 0; i < length; ++i){ 
                if (letter == 0) {
                    if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)]) {
                        ret.matches += 1;
                    }
                    ++q;
                    ++p;
                } else {
                    if (letter == 1) ++p;
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

}; /* namespace pgraph */

