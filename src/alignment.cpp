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
#include "alignment_defs.hpp"
#include "blosum/blosum40.h"
#include "blosum/blosum45.h"
#include "blosum/blosum50.h"
#include "blosum/blosum62.h"
#include "blosum/blosum75.h"
#include "blosum/blosum80.h"
#include "blosum/blosum90.h"
#include "ssw.h"

#define CROW(i) ((i)%2)
#define PROW(i) ((i-1)%2)
#define BLOSUM(ch1, ch2) blosum[MAP_BLOSUM[(ch1)-'A']][MAP_BLOSUM[(ch2)-'A']]
#define SCORE(sub, map, first, ch1, ch2) \
    (sub)[(map)[(ch1)-(first)]][(map)[(ch2)-(first)]]
#define NEG_ADD(x, y) (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y))
#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))

using std::cerr;
using std::endl;


namespace pgraph {


/** points to selected blosum table (default blosum62) */
static const int (* restrict blosum)[24] = blosum62;
static const int8_t * restrict blosum__ = blosum62__;


cell_t **allocate_cell_table(size_t nrow, size_t ncol)
{
    size_t i = 0;
    cell_t **table = NULL;

    table = new cell_t*[nrow];
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        table[i] = new cell_t[ncol + 1];
    }

    return table;
}


int **allocate_int_table(size_t nrow, size_t ncol)
{
    size_t i = 0;
    int **table = NULL;

    table = new int*[nrow];
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        table[i] = new int[ncol + 1];
    }

    return table;
}


void free_cell_table(cell_t **table, size_t nrow)
{
    size_t i;
    for (i = 0; i < nrow; i++) {
        delete [] table[i];
    }
    delete [] table;
}


void free_int_table(int **table, size_t nrow)
{
    size_t i;
    for (i = 0; i < nrow; i++) {
        delete [] table[i];
    }
    delete [] table;
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


int self_score(const char * const restrict seq, size_t len, int match)
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
            blosum__ = blosum40__;
            break;
        case 45:
            blosum = blosum45;
            blosum__ = blosum45__;
            break;
        case 62:
            blosum = blosum62;
            blosum__ = blosum62__;
            break;
        case 75:
            blosum = blosum75;
            blosum__ = blosum75__;
            break;
        case 80:
            blosum = blosum80;
            blosum__ = blosum80__;
            break;
        case 90:
            blosum = blosum90;
            blosum__ = blosum90__;
            break;
        default:
            cerr << "invalid blosum number (" << number << ")" << endl;
            assert(0);
    }
}


#define AFFINE_GAP_ALIGN_DECL                                               \
    size_t i = 0;                       /* first sequence char index */     \
    size_t j = 0;                       /* second sequence char index */    \
    cell_t lastCol = {INT_MIN, 0, 0};   /* max cell in last col */          \
    cell_t lastRow = {INT_MIN, 0, 0};   /* max cell in last row */          \
    cell_t * restrict tI = NULL;        /* current DP table row */          \
    cell_t * restrict pI = NULL;        /* previous DP table row */

#define AFFINE_GAP_ALIGN_ASSERT \
    assert(s1);                 \
    assert(s2);                 \
    assert(s1Len > 0);          \
    assert(s2Len > 0);          \
    assert(tbl);                \
    assert(del);                \
    assert(ins);

#define AFFINE_GAP_ALIGN_BODY(CALCULATE_SCORE)                              \
    /* init first row of 3 tables */                                        \
    tbl[0][0].score = 0;                                                    \
    tbl[0][0].matches = 0;                                                  \
    tbl[0][0].length = 0;                                                   \
    del[0][0] = 0;                                                          \
    ins[0][0] = 0;                                                          \
    tI = tbl[0];                                                            \
    for (j = 1; j <= s2Len; j++) {                                          \
        tI[j].score = open + j * gap;                                       \
        tI[j].matches = 0;                                                  \
        tI[j].length = 0;                                                   \
        del[0][j] = INT_MIN;                                                \
        ins[0][j] = open + j * gap;                                         \
    }                                                                       \
                                                                            \
    for (i = 1; i <= s1Len; ++i) {                                          \
        int dig = 0;        /* diagonal value in DP table */                \
        int up = 0;         /* upper value in DP table */                   \
        int left = 0;       /* left value in DP table */                    \
        int maxScore = 0;   /* max of dig, up, and left in DP table */      \
        size_t cr = 0;      /* current row index in DP table */             \
        size_t pr = 0;      /* previous row index in DP table */            \
        char ch1 = 0;       /* current character in first sequence */       \
                                                                            \
        cr = CROW(i);                                                       \
        pr = PROW(i);                                                       \
        ch1 = s1[i - 1];                                                    \
        tI = tbl[cr];                                                       \
        pI = tbl[pr];                                                       \
                                                                            \
        /* init first column of 3 tables */                                 \
        tI[0].score = open + i * gap;                                       \
        tI[0].matches = 0;                                                  \
        tI[0].length = 0;                                                   \
        del[cr][0] = open + i * gap;                                        \
        ins[cr][0] = INT_MIN;                                               \
                                                                            \
        for (j = 1; j <= s2Len; j++) {                                      \
            int tmp1 = 0;   /* temporary during DP calculation */           \
            int tmp2 = 0;   /* temporary during DP calculation */           \
            char ch2 = 0;   /* current character in second sequence */      \
                                                                            \
            ch2 = s2[j - 1];                                                \
                                                                            \
            /* overflow could happen, INT_MIN-1 = 2147483647 */             \
            tmp1 = pI[j].score + open + gap;                                \
            tmp2 = NEG_ADD(del[pr][j], gap);                                \
            up = MAX(tmp1, tmp2);                                           \
            del[cr][j] = up;                                                \
            tmp1 = tI[j - 1].score + open + gap;                            \
            tmp2 = NEG_ADD(ins[cr][j - 1], gap);                            \
            left = MAX(tmp1, tmp2);                                         \
            ins[cr][j] = left;                                              \
            maxScore = MAX(up, left);                                       \
            dig = pI[j - 1].score + CALCULATE_SCORE;                        \
            maxScore = MAX(maxScore, dig);                                  \
            tI[j].score = maxScore;                                         \
                                                                            \
            if (maxScore == dig) {                                          \
                tI[j].matches = pI[j - 1].matches + ((ch1 == ch2) ? 1 : 0); \
                tI[j].length = pI[j - 1].length + 1;                        \
            }                                                               \
            else if (maxScore == up) {                                      \
                tI[j].matches = pI[j].matches;                              \
                tI[j].length = pI[j].length + 1;                            \
            }                                                               \
            else {                                                          \
                tI[j].matches = tI[j - 1].matches;                          \
                tI[j].length = tI[j - 1].length + 1;                        \
            }                                                               \
                                                                            \
            /* track the maximum of last row */                             \
            if (i == s1Len && tI[j].score > lastRow.score) {                \
                lastRow = tI[j];                                            \
            }                                                               \
        } /* end of j loop */                                               \
                                                                            \
        assert(j == (s2Len + 1));                                           \
                                                                            \
        /* update the maximum of last column */                             \
        if (tI[s2Len].score > lastCol.score) {                              \
            lastCol = tI[s2Len];                                            \
        }                                                                   \
    } /* end of i loop */                                                   \
                                                                            \
    return (lastCol.score > lastRow.score) ? lastCol : lastRow;


cell_t affine_gap_align(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first)
{
    AFFINE_GAP_ALIGN_DECL

    AFFINE_GAP_ALIGN_ASSERT
    assert(sub);
    assert(map);

    AFFINE_GAP_ALIGN_BODY(SCORE(sub, map, first, ch1, ch2))
}


cell_t affine_gap_align(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        match_t callback)
{
    AFFINE_GAP_ALIGN_DECL

    AFFINE_GAP_ALIGN_ASSERT

    AFFINE_GAP_ALIGN_BODY(callback(ch1, ch2))
}


cell_t affine_gap_align(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        int match,
        int mismatch)
{
    AFFINE_GAP_ALIGN_DECL

    AFFINE_GAP_ALIGN_ASSERT

    AFFINE_GAP_ALIGN_BODY((ch1 == ch2 ? match : mismatch))
}


cell_t affine_gap_align_blosum(
        const char * const restrict s1, size_t s1Len,
        const char * const restrict s2, size_t s2Len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap)
{
    AFFINE_GAP_ALIGN_DECL

    AFFINE_GAP_ALIGN_ASSERT

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
        tI[0].score = open + j * gap;
        tI[0].matches = 0;                                                  
        tI[0].length = 0;                                                   
        del[cr][0] = open + i * gap;                                        
        ins[cr][0] = INT_MIN;                                               
                                                                            
        for (j = 1; j <= s2Len; j++) {                                      
            int tmp1 = 0;   /* temporary during DP calculation */           
            int tmp2 = 0;   /* temporary during DP calculation */           
            int tmp3 = 0;   /* temporary during DP calculation */           
            char ch2 = 0;   /* current character in second sequence */      
                                                                            
            ch2 = s2[j - 1];                                                
                                                                            
            /* overflow could happen, INT_MIN-1 = 2147483647 */             
            tmp1 = NEG_ADD(NEG_ADD(pI[j].score, open), gap);
            tmp2 = NEG_ADD(del[pr][j], gap);                                
            up = MAX(tmp1, tmp2);                                           
            del[cr][j] = up;                                                
            tmp1 = NEG_ADD(NEG_ADD(tI[j - 1].score, open), gap);
            tmp2 = NEG_ADD(ins[cr][j - 1], gap);                            
            left = MAX(tmp1, tmp2);                                         
            ins[cr][j] = left;                                              
            maxScore = MAX(up, left);
            tmp3 = BLOSUM(ch1, ch2);                       
            dig = NEG_ADD(pI[j - 1].score, tmp3);
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
                                                                            
    return (lastCol.score > lastRow.score) ? lastCol : lastRow;
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


bool is_edge_blosum(
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


cell_t affine_gap_align_blosum_ssw(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open, int gap)
{
    cell_t ret;
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
    result = ssw_align(profile, s2_num, s2_len, open, gap, 2, 0, 0, s1_len/2);

    ret.score = result->score1;
    ret.matches = 0;
    ret.length = 0;

    s_align* a = result;
    const char * const restrict ref_seq = s2;
    const char * const restrict read_seq = s1;
    const int8_t * const restrict mat = blosum__;
    int32_t n = 24;

    if (a->cigar) {
        int32_t i, c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
        while (e < a->cigarLen || left > 0) {
            int32_t count = 0;
            int32_t q = qb;
            int32_t p = pb;
            for (c = e; c < a->cigarLen; ++c) {
                int32_t letter = 0xf&*(a->cigar + c);
                int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                int32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) {
                    if (letter == 1) {}
                    else {
                        ++ q;
                    }
                    ++ count;
                    if (count == 60) goto step2;
                }
            }
step2:
            q = qb;
            count = 0;
            for (c = e; c < a->cigarLen; ++c) {
                int32_t letter = 0xf&*(a->cigar + c);
                int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                int32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i){ 
                    if (letter == 0) {
                        if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)]){
                            // |
                            ret.matches += 1;
                            ret.length += 1;
                        }
                        else if (mat[table[(int)*(ref_seq + q)]*n+table[(int)*(read_seq + p)]]>0){
                            // :
                            ret.length += 1;
                        }
                        else {
                            // *
                            ret.length += 1;
                        }
                        ++q;
                        ++p;
                    } else {
                        if (letter == 1) ++p;
                        else ++q;
                    }
                    ++ count;
                    if (count == 60) {
                        qb = q;
                        goto step3;
                    }
                }
            }
step3:
            p = pb;
            count = 0;
            for (c = e; c < a->cigarLen; ++c) {
                int32_t letter = 0xf&*(a->cigar + c);
                int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                int32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) { 
                    if (letter == 2) {}
                    else {
                        ++p;
                    }
                    ++ count;
                    if (count == 60) {
                        pb = p;
                        left = l - i - 1;
                        e = (left == 0) ? (c + 1) : c;
                        goto end;
                    }
                }
            }
            e = c;
            left = 0;
end:
            ;
        }
    }

    init_destroy(profile);
    if (result->cigar) free(result->cigar);
    free(result);
    delete [] s1_num;
    delete [] s2_num;

    return ret;
}

}; /* namespace pgraph */

