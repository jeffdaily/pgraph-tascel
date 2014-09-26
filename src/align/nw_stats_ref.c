#include "config.h"

#include <stdlib.h>

#ifdef ALIGN_EXTRA
#include "align/align_debug.h"
#else
#include "align/align.h"
#endif
#include "blosum/blosum_map.h"

#ifdef ALIGN_EXTRA
#define FNAME nw_stats_debug
#else
#define FNAME nw_stats
#endif

int FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * matches, int * length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
        , int * const restrict match_table
        , int * const restrict length_table
#endif
        )
{
    int * const restrict s1 = (int * const restrict)malloc(sizeof(int)*s1Len);
    int * const restrict s2 = (int * const restrict)malloc(sizeof(int)*s2Len);
    int i = 0;
    int j = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }

    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF_32;
    mch_pr[0] = 0;
    len_pr[0] = 0;
    
    /* first row */
    for (j=1; j<=s2Len; ++j) {
        tbl_pr[j] = -open -(j-1)*gap;
        del_pr[j] = NEG_INF_32;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    /* iter over first sequence */
    for (i=1; i<=s1Len; ++i) {
        const int * const restrict matrow = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nlength = len_pr[0];
        int Wscore = -open - (i-1)*gap;
        int Wmatches = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        len_pr[0] = Wlength;
        for (j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + matrow[s2[j-1]];
            if ((tbl_pr[j] >= del_pr[j]) && (tbl_pr[j] >= ins_cr)) {
                Wscore = tbl_pr[j];
                Wmatches  = NWmatches + (s1[i-1] == s2[j-1]);
                Wlength  = NWlength + 1;
            } else if (del_pr[j] >= ins_cr) {
                Wscore = del_pr[j];
                Wmatches  = Nmatches;
                Wlength  = Nlength + 1;
            } else {
                Wscore = ins_cr;
                Wmatches  = Wmatches;
                Wlength  = Wlength + 1;
            }
            tbl_pr[j] = Wscore;
            mch_pr[j] = Wmatches;
            len_pr[j] = Wlength;
#ifdef ALIGN_EXTRA
            score_table[(i-1)*s2Len + (j-1)] = Wscore;
            match_table[(i-1)*s2Len + (j-1)] = Wmatches;
            length_table[(i-1)*s2Len + (j-1)] = Wlength;
#endif
        }
    }

    free(s1);
    free(s2);
    *matches = mch_pr[s2Len];
    *length = len_pr[s2Len];
    return tbl_pr[s2Len];
}

