#include "config.h"

#include <cassert>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>

#include <emmintrin.h>

#include "blosum/blosum62.h"
#include "timer.h"

using namespace ::std;

#define DEBUG 0
#define DEBUG_MATCHES 0
#define DEBUG_LENGTH 0
#define SHORT_TEST 0

/* This table is used to transform amino acid letters into numbers. */
static const int MAP_BLOSUM_[256] = {
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 22, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 22, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23
};

#define BLOSUM(ch1, ch2) (matrix[MAP_BLOSUM_[(ch1)]][MAP_BLOSUM_[(ch2)]])
#define BLOSUM_(ch1, ch2) (matrix[(ch1)][(ch2)])

#define BLOSUM0_(ch1, ch2) (matrow0[(ch2)])
#define BLOSUM1_(ch1, ch2) (matrow1[(ch2)])
#define BLOSUM2_(ch1, ch2) (matrow2[(ch2)])
#define BLOSUM3_(ch1, ch2) (matrow3[(ch2)])
#define BLOSUM4_(ch1, ch2) (matrow4[(ch2)])
#define BLOSUM5_(ch1, ch2) (matrow5[(ch2)])
#define BLOSUM6_(ch1, ch2) (matrow6[(ch2)])
#define BLOSUM7_(ch1, ch2) (matrow7[(ch2)])

#define MAX(a,b) ((a)>(b)?(a):(b))

/* on OSX, _mm_extract_epi16 got wrong answer, but taking the union of
 * the vector and extracting that way seemed to work... */
#define EXTRACT extract
//#define EXTRACT _mm_extract_epi16

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
                    (m) = EXTRACT((vm), 0)

struct DP_t {
    int score;
    int matches;
    int length;
    DP_t() : score(0), matches(0), length(0) {}
    DP_t(int score, int matches, int length)
        : score(score), matches(matches), length(length) {}
};

//static const int16_t NEG_INF = SHRT_MIN / 2;
static const int16_t NEG_INF = -50;

static inline int extract(const __m128i &m, const int &pos)
{
    union {
        __m128i m;
        int16_t v[8];
    } tmp;
    tmp.m = m;
    return tmp.v[pos];
}


static void init_vect(int len, int *vec, int val)
{
    int i;
    for (i=0; i<len; ++i) {
        vec[i] = val;
    }
}


#if DEBUG
static void print_array(
        const int * const restrict array,
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int pad)
{
    int i;
    int j;
    int padb = pad/2;

    printf(" ");
    for (j=0; j<pad/2; ++j) {
        printf("   $");
    }
    for (j=0; j<lena; ++j) {
        printf("   %c", seqA[j]);
    }
    for (j=0; j<pad/2; ++j) {
        printf("   $");
    }
    printf("\n");
    for (i=0; i<lenb; ++i) {
        printf("%c", seqB[i]);
        for (j=0; j<lena+pad; ++j) {
            printf(" %3d", array[i*(lena+pad) + j]);
        }
        printf("\n");
    }
    for (i=lenb; i<(lenb+padb); ++i) {
        printf("$");
        for (j=0; j<lena+pad; ++j) {
            printf(" %3d", array[i*(lena+pad) + j]);
        }
        printf("\n");
    }
}


static void print_array2(
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int pad)
{
    int i;
    int j;
    int padb = pad/2;

    printf(" ");
    for (j=0; j<pad/2; ++j) {
        printf("   $");
    }
    printf("   *");
    for (j=0; j<s2Len; ++j) {
        printf("   %c", s2[j]);
    }
    for (j=0; j<pad/2; ++j) {
        printf("   $");
    }
    printf("\n");
    {
        printf("*");
        for (j=0; j<s2Len+pad+1; ++j) {
            printf(" %3d", array[j]);
        }
        printf("\n");
    }
    for (i=0; i<s1Len; ++i) {
        printf("%c", s1[i]);
        for (j=0; j<s2Len+pad+1; ++j) {
            printf(" %3d", array[(i+1)*(s2Len+pad+1) + j]);
        }
        printf("\n");
    }
    for (i=s1Len; i<(s1Len+padb); ++i) {
        printf("$");
        for (j=0; j<s2Len+pad+1; ++j) {
            printf(" %3d", array[(i+1)*(s2Len+pad+1) + j]);
        }
        printf("\n");
    }
}


static void print_array2(
        const int * const restrict array,
        const int * const restrict s1, const int s1Len,
        const int * const restrict s2, const int s2Len,
        const int pad)
{
    int i;
    int j;
    int padb = pad/2;

    printf("  ");
    for (j=0; j<pad/2; ++j) {
        printf("   $");
    }
    printf("   *");
    for (j=0; j<s2Len; ++j) {
        printf("  %2d", s2[j]);
    }
    for (j=0; j<pad/2; ++j) {
        printf("   $");
    }
    printf("\n");
    {
        printf(" *");
        for (j=0; j<s2Len+pad+1; ++j) {
            printf(" %3d", array[j]);
        }
        printf("\n");
    }
    for (i=0; i<s1Len; ++i) {
        printf("%2d", s1[i]);
        for (j=0; j<s2Len+pad+1; ++j) {
            printf(" %3d", array[(i+1)*(s2Len+pad+1) + j]);
        }
        printf("\n");
    }
    for (i=s1Len; i<(s1Len+padb); ++i) {
        printf(" $");
        for (j=0; j<s2Len+pad+1; ++j) {
            printf(" %3d", array[(i+1)*(s2Len+pad+1) + j]);
        }
        printf("\n");
    }
}


static void print_m128i_16(const char *name, const __m128i &m) {
    union {
        __m128i m;
        int16_t v[8];
    } tmp;
    tmp.m = m;
    printf("%s={%3d,%3d,%3d,%3d,%3d,%3d,%3d,%3d}\n", name,
            int(tmp.v[0]), int(tmp.v[1]), int(tmp.v[2]), int(tmp.v[3]),
            int(tmp.v[4]), int(tmp.v[5]), int(tmp.v[6]), int(tmp.v[7]));
}
#endif


static int sw_woz(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i = 0; 
    int j = 0;
    int last_nogap = 0;
    int prev_nogap = 0;
    int score = 0;
    init_vect(lena, nogap, 0);
    init_vect(lena, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+0)*lenb(=%d) = %d\n",
            lena, lenb, (lena+0)*lenb);
    int *array = new int[(lena+0)*lenb];
    init_vect((lena+0)*lenb, array, -9);
#endif
    for (i=0; i<lenb; ++i) {
        int a_gap;
        last_nogap = prev_nogap = 0;
        a_gap = -gap_open;
        for (j=0; j<lena; ++j) {
            a_gap = MAX((last_nogap - gap_open), (a_gap - gap_ext));
            b_gap[j] = MAX((nogap[j] - gap_open), (b_gap[j] - gap_ext));
            last_nogap = MAX((prev_nogap + BLOSUM(seqA[j],seqB[i])), 0);
            last_nogap = MAX(last_nogap, a_gap);
            last_nogap = MAX(last_nogap, b_gap[j]);
            prev_nogap = nogap[j];
            nogap[j] = last_nogap;
#if DEBUG
            array[i*lena + j] = last_nogap;
#endif
            score = MAX(score, last_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 0);
    delete [] array;
#endif
    return score;
}


static int nw(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF;
#if DEBUG
    printf(" %3d", tbl_pr[0]);
#endif
    
    /* first row */
    for (int j=1; j<=s2Len; ++j) {
        tbl_pr[j] = -open -(j-1)*gap;
        del_pr[j] = NEG_INF;
#if DEBUG
        printf(" %3d", tbl_pr[j]);
#endif
    }
#if DEBUG
    printf("\n");
#endif

    /* iter over first sequence */
    for (int i=1; i<=s1Len; ++i) {
        const int * const restrict matrow = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = -open - (i-1)*gap;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
#if DEBUG
        printf(" %3d", Wscore);
#endif
        for (int j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            //tbl_pr[j] = NWscore + BLOSUM_(s1[i-1],s2[j-1]);
            tbl_pr[j] = NWscore + matrow[s2[j-1]];
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
#if DEBUG
            printf(" %3d", Wscore);
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }

    delete [] s1;
    delete [] s2;
    return tbl_pr[s2Len];
}


static DP_t nw_stats(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF;
    mch_pr[0] = 0;
    len_pr[0] = 0;
#if DEBUG
#if DEBUG_MATCHES
    printf(" %3d", mch_pr[0]);
#elif DEBUG_LENGTH
    printf(" %3d", len_pr[0]);
#else
    printf(" %3d", tbl_pr[0]);
#endif
#endif
    
    /* first row */
    for (int j=1; j<=s2Len; ++j) {
        tbl_pr[j] = -open -(j-1)*gap;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        printf(" %3d", mch_pr[j]);
#elif DEBUG_LENGTH
        printf(" %3d", len_pr[j]);
#else
        printf(" %3d", tbl_pr[j]);
#endif
#endif
    }
#if DEBUG
    printf("\n");
#endif

    /* iter over first sequence */
    for (int i=1; i<=s1Len; ++i) {
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nlength = len_pr[0];
        int Wscore = -open - (i-1)*gap;
        int Wmatches = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        len_pr[0] = Wlength;
#if DEBUG
#if DEBUG_MATCHES
        printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
        printf(" %3d", Wlength);
#else
        printf(" %3d", Wscore);
#endif
#endif
        for (int j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM_(s1[i-1],s2[j-1]);
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
#if DEBUG
#if DEBUG_MATCHES
            printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
            printf(" %3d", Wlength);
#else
            printf(" %3d", Wscore);
#endif
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }

    delete [] s1;
    delete [] s2;
    return DP_t(tbl_pr[s2Len], mch_pr[s2Len], len_pr[s2Len]);
}


static int sg(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int score = NEG_INF;
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF;
#if DEBUG
    printf(" %3d", tbl_pr[0]);
#endif
    
    /* first row */
    for (int j=1; j<=s2Len; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
#if DEBUG
        printf(" %3d", tbl_pr[j]);
#endif
    }
#if DEBUG
    printf("\n");
#endif

    /* iter over first sequence */
    for (int i=1; i<s1Len; ++i) {
        const int * const restrict matrow0 = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = 0;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
#if DEBUG
        printf(" %3d", Wscore);
#endif
        for (int j=1; j<s2Len; ++j) {
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM0_(s1[i-1],s2[j-1]);
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
#if DEBUG
            printf(" %3d", Wscore);
#endif
        }
        {
            int j=s2Len;
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM0_(s1[i-1],s2[j-1]);
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
            score = MAX(score, Wscore);
#if DEBUG
            printf(" %3d", Wscore);
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }
    {
        int i=s1Len;
        const int * const restrict matrow0 = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = 0;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
#if DEBUG
        printf(" %3d", Wscore);
#endif
        for (int j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM0_(s1[i-1],s2[j-1]);
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
            score = MAX(score, Wscore);
#if DEBUG
            printf(" %3d", Wscore);
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }

    delete [] s1;
    delete [] s2;
    return score;
}


static DP_t sg_stats(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    int score = NEG_INF;
    int matches = NEG_INF;
    int length = NEG_INF;
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF;
    mch_pr[0] = 0;
    len_pr[0] = 0;
#if DEBUG
#if DEBUG_MATCHES
    printf(" %3d", mch_pr[0]);
#elif DEBUG_LENGTH
    printf(" %3d", len_pr[0]);
#else
    printf(" %3d", tbl_pr[0]);
#endif
#endif
    
    /* first row */
    for (int j=1; j<=s2Len; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        printf(" %3d", mch_pr[j]);
#elif DEBUG_LENGTH
        printf(" %3d", len_pr[j]);
#else
        printf(" %3d", tbl_pr[j]);
#endif
#endif
    }
#if DEBUG
    printf("\n");
#endif

    /* iter over first sequence */
    for (int i=1; i<s1Len; ++i) {
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nlength = len_pr[0];
        int Wscore = 0;
        int Wmatches = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        len_pr[0] = Wlength;
#if DEBUG
#if DEBUG_MATCHES
        printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
        printf(" %3d", Wlength);
#else
        printf(" %3d", Wscore);
#endif
#endif
        for (int j=1; j<s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM_(s1[i-1],s2[j-1]);
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
#if DEBUG
#if DEBUG_MATCHES
            printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
            printf(" %3d", Wlength);
#else
            printf(" %3d", Wscore);
#endif
#endif
        }
        {
            int j=s2Len;
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM_(s1[i-1],s2[j-1]);
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
            if (Wscore > score) {
                score = Wscore;
                matches = Wmatches;
                length = Wlength;
            }
#if DEBUG
#if DEBUG_MATCHES
            printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
            printf(" %3d", Wlength);
#else
            printf(" %3d", Wscore);
#endif
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }
    {
        int i=s1Len;
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nlength = len_pr[0];
        int Wscore = 0;
        int Wmatches = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        len_pr[0] = Wlength;
#if DEBUG
#if DEBUG_MATCHES
        printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
        printf(" %3d", Wlength);
#else
        printf(" %3d", Wscore);
#endif
#endif
        for (int j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM_(s1[i-1],s2[j-1]);
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
            if (Wscore > score) {
                score = Wscore;
                matches = Wmatches;
                length = Wlength;
            }
#if DEBUG
#if DEBUG_MATCHES
            printf(" %3d", Wmatches);
#elif DEBUG_LENGTH
            printf(" %3d", Wlength);
#else
            printf(" %3d", Wscore);
#endif
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }

    delete [] s1;
    delete [] s2;
    return DP_t(score, matches, length);
}


static int sw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int score = NEG_INF;
    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF;
#if DEBUG
    printf(" %3d", tbl_pr[0]);
#endif
    
    /* first row */
    for (int j=1; j<=s2Len; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
#if DEBUG
        printf(" %3d", tbl_pr[j]);
#endif
    }
#if DEBUG
    printf("\n");
#endif

    /* iter over first sequence */
    for (int i=1; i<=s1Len; ++i) {
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = 0;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
#if DEBUG
        printf(" %3d", Wscore);
#endif
        for (int j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM(s1[i-1],s2[j-1]);
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
            score = MAX(score,Wscore);
#if DEBUG
            printf(" %3d", Wscore);
#endif
        }
#if DEBUG
        printf("\n");
#endif
    }

    return score;
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


static int nw_sse8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int score = 0;
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
#if DEBUG
    printf("array length (s2Len(=%d)+14+1)*(s1Len(=%d)+7+1) = %d\n",
            s2Len, s1Len, (s2Len+14+1)*(s1Len+7+1));
    int *array = new int[(s2Len+14+1)*(s1Len+7+1)];
    init_vect((s2Len+14+1)*(s1Len+7+1), array, -99);
#endif

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
#if DEBUG
    array[7] = tbl_pr[7];
#endif
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = -open -(j-8)*gap;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (int i=1; i<=s1Len; i+=8) {
        int j;
        __m128i NWscore = _mm_set1_epi16(NEG_INF);
        __m128i Nscore  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl    = _mm_set1_epi16(NEG_INF);
        __m128i vDel    = _mm_set1_epi16(NEG_INF);
        __m128i vIns    = _mm_set1_epi16(NEG_INF);
        __m128i vMat;

        const int * const restrict matrow0 = matrix[s1[i-1+0]];
        const int * const restrict matrow1 = matrix[s1[i-1+1]];
        const int * const restrict matrow2 = matrix[s1[i-1+2]];
        const int * const restrict matrow3 = matrix[s1[i-1+3]];
        const int * const restrict matrow4 = matrix[s1[i-1+4]];
        const int * const restrict matrow5 = matrix[s1[i-1+5]];
        const int * const restrict matrow6 = matrix[s1[i-1+6]];
        const int * const restrict matrow7 = matrix[s1[i-1+7]];

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, -open-(i-1)*gap, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = extract(Wscore, 7);
#endif

        /* j = 1 */
        j = 1;
        NWscore = Nscore;
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);
        vDel    = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+0)*gap, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = extract(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = extract(Wscore, 6);
#endif

        /* j = 2 */
        j = 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+1)*gap, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif

        /* j = 3 */
        j = 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+2)*gap, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif

        /* j = 4 */
        j = 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+3)*gap, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif

        /* j = 5 */
        j = 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+4)*gap, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif

        /* j = 6 */
        j = 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+5)*gap, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wscore, 1);
#endif

        /* j = 7 */
        j = 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+6)*gap, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = -open-(i+6)*gap;
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wscore, 0);
#endif

        for (j=8; j<=s2Len; ++j) {
            NWscore = Nscore;
            Nscore = vshift16(Wscore, tbl_pr[j+7]);
            vDel = vshift16(vDel, del_pr[j+7]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(Nscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(Wscore,vOpen),
                    _mm_sub_epi16(vIns,vGap));
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
            vTbl = _mm_add_epi16(NWscore, vMat);
            vTbl = _mm_max_epi16(vTbl, vDel);
            vTbl = _mm_max_epi16(vTbl, vIns);
            Wscore = vTbl;
            tbl_pr[j] = EXTRACT(vTbl, 0);
            del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        }
        if (i+0 == s1Len) score = EXTRACT(vTbl, 7);

        /* j = s2Len + 1 */
        j = s2Len + 1;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+1 == s1Len) score = EXTRACT(vTbl, 6);

        /* j = s2Len + 2 */
        j = s2Len + 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+2 == s1Len) score = EXTRACT(vTbl, 5);

        /* j = s2Len + 3 */
        j = s2Len + 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+3 == s1Len) score = EXTRACT(vTbl, 4);

        /* j = s2Len + 4 */
        j = s2Len + 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+4 == s1Len) score = EXTRACT(vTbl, 3);

        /* j = s2Len + 5 */
        j = s2Len + 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+5 == s1Len) score = EXTRACT(vTbl, 2);

        /* j = s2Len + 6 */
        j = s2Len + 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+6 == s1Len) score = EXTRACT(vTbl, 1);

        /* j = s2Len + 7 */
        j = s2Len + 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+7 == s1Len) score = EXTRACT(vTbl, 0);
    }

#if DEBUG
    print_array2(array, s1, s1Len, s2, s2Len, 14);
    delete [] array;
#endif
    delete [] s1;
    delete [] s2;
    return score;
}


static DP_t nw_stats_sse8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    int score = 0;
    int match = 0;
    int length = 0;
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
#if DEBUG
    printf("array length (s2Len(=%d)+14+1)*(s1Len(=%d)+7+1) = %d\n",
            s2Len, s1Len, (s2Len+14+1)*(s1Len+7+1));
    int *array = new int[(s2Len+14+1)*(s1Len+7+1)];
    init_vect((s2Len+14+1)*(s1Len+7+1), array, -99);
#endif

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        array[j] = mch_pr[j];
#elif DEBUG_LENGTH
        array[j] = len_pr[j];
#else
        array[j] = tbl_pr[j];
#endif
#endif
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
    mch_pr[7] = 0;
    len_pr[7] = 0;
#if DEBUG
#if DEBUG_MATCHES
    array[7] = mch_pr[7];
#elif DEBUG_LENGTH
    array[7] = len_pr[7];
#else
    array[7] = tbl_pr[7];
#endif
#endif
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = -open -(j-8)*gap;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        array[j] = mch_pr[j];
#elif DEBUG_LENGTH
        array[j] = len_pr[j];
#else
        array[j] = tbl_pr[j];
#endif
#endif
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        array[j] = mch_pr[j];
#elif DEBUG_LENGTH
        array[j] = len_pr[j];
#else
        array[j] = tbl_pr[j];
#endif
#endif
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

        const int * const restrict matrow0 = matrix[s1[i-1+0]];
        const int * const restrict matrow1 = matrix[s1[i-1+1]];
        const int * const restrict matrow2 = matrix[s1[i-1+2]];
        const int * const restrict matrow3 = matrix[s1[i-1+3]];
        const int * const restrict matrow4 = matrix[s1[i-1+4]];
        const int * const restrict matrow5 = matrix[s1[i-1+5]];
        const int * const restrict matrow6 = matrix[s1[i-1+6]];
        const int * const restrict matrow7 = matrix[s1[i-1+7]];

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength,7);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wscore, 7);
#endif
#endif

        /* j = 1 */
        j = 1;
        /* this block never changes, so let's not copy-and-paste */
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wscore, 6);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wscore, 1);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wscore, 0);
#endif
#endif

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
#if DEBUG
#if DEBUG_MATCHES
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
        if (i+7 == s1Len) {
            score = EXTRACT(vTbl, 0);
            match = EXTRACT(Wmatch, 0);
            length= EXTRACT(Wlength, 0);
        }
    }

#if DEBUG
    print_array2(array, s1, s1Len, s2, s2Len, 14);
    delete [] array;
#endif
    delete [] s1;
    delete [] s2;
    return DP_t(score, match, length);
}


/* global alignment without end gap penalties */
static int sg_sse8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int score = NEG_INF;
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
#if DEBUG
    printf("array length (s2Len(=%d)+14+1)*(s1Len(=%d)+7+1) = %d\n",
            s2Len, s1Len, (s2Len+14+1)*(s1Len+7+1));
    int *array = new int[(s2Len+14+1)*(s1Len+7+1)];
    init_vect((s2Len+14+1)*(s1Len+7+1), array, -99);
#endif

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
#if DEBUG
    array[7] = tbl_pr[7];
#endif
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (int i=1; i<=s1Len; i+=8) {
        bool last_pass = (i+8>=s1Len);
        int offset = 7 - (s1Len - i);
        int j;
        __m128i NWscore = _mm_set1_epi16(NEG_INF);
        __m128i Nscore  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl    = _mm_set1_epi16(NEG_INF);
        __m128i vDel    = _mm_set1_epi16(NEG_INF);
        __m128i vIns    = _mm_set1_epi16(NEG_INF);
        __m128i vMat;

        const int * const restrict matrow0 = matrix[s1[i-1+0]];
        const int * const restrict matrow1 = matrix[s1[i-1+1]];
        const int * const restrict matrow2 = matrix[s1[i-1+2]];
        const int * const restrict matrow3 = matrix[s1[i-1+3]];
        const int * const restrict matrow4 = matrix[s1[i-1+4]];
        const int * const restrict matrow5 = matrix[s1[i-1+5]];
        const int * const restrict matrow6 = matrix[s1[i-1+6]];
        const int * const restrict matrow7 = matrix[s1[i-1+7]];

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, 0, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wscore, 7);
#endif

        /* j = 1 */
        j = 1;
        NWscore = Nscore;
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);
        vDel    = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wscore, 6);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = 2 */
        j = 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = 3 */
        j = 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = 4 */
        j = 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = 5 */
        j = 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = 6 */
        j = 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wscore, 1);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = 7 */
        j = 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = 0;
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wscore, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        if (j == s2Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        for (j=8; j<=s2Len; ++j) {
            NWscore = Nscore;
            Nscore = vshift16(Wscore, tbl_pr[j+7]);
            vDel = vshift16(vDel, del_pr[j+7]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(Nscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(Wscore,vOpen),
                    _mm_sub_epi16(vIns,vGap));
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
            vTbl = _mm_add_epi16(NWscore, vMat);
            vTbl = _mm_max_epi16(vTbl, vDel);
            vTbl = _mm_max_epi16(vTbl, vIns);
            Wscore = vTbl;
            tbl_pr[j] = EXTRACT(vTbl, 0);
            del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
            if (last_pass) {
                int tmp = vextract16(vTbl, offset);
                score = MAX(score,tmp);
            }
            if (j == s2Len) {
                int tmp = EXTRACT(vTbl, 7);
                score = MAX(score,tmp);
            }
        }

        /* j = s2Len + 1 */
        j = s2Len + 1;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 6);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 2 */
        j = s2Len + 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 5);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 3 */
        j = s2Len + 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 4);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 4 */
        j = s2Len + 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 3);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 5 */
        j = s2Len + 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 2);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 6 */
        j = s2Len + 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 1);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 7 */
        j = s2Len + 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
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
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (last_pass) {
            int tmp = vextract16(vTbl, offset);
            score = MAX(score,tmp);
        }
        {
            int tmp = EXTRACT(vTbl, 0);
            score = MAX(score,tmp);
        }
    }

#if DEBUG
    print_array2(array, s1, s1Len, s2, s2Len, 14);
    delete [] array;
#endif
    delete [] s1;
    delete [] s2;
    return score;
}


static DP_t sg_stats_sse8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr)
{
    int score = NEG_INF;
    int match = 0;
    int length = 0;
    int * const restrict s1 = new int[s1Len];
    int * const restrict s2 = new int[s2Len];
    for (int i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (int j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }
#if DEBUG
    printf("array length (s2Len(=%d)+14+1)*(s1Len(=%d)+7+1) = %d\n",
            s2Len, s1Len, (s2Len+14+1)*(s1Len+7+1));
    int *array = new int[(s2Len+14+1)*(s1Len+7+1)];
    init_vect((s2Len+14+1)*(s1Len+7+1), array, -99);
#endif

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        array[j] = mch_pr[j];
#elif DEBUG_LENGTH
        array[j] = len_pr[j];
#else
        array[j] = tbl_pr[j];
#endif
#endif
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
    mch_pr[7] = 0;
    len_pr[7] = 0;
#if DEBUG
#if DEBUG_MATCHES
    array[7] = mch_pr[7];
#elif DEBUG_LENGTH
    array[7] = len_pr[7];
#else
    array[7] = tbl_pr[7];
#endif
#endif
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        array[j] = mch_pr[j];
#elif DEBUG_LENGTH
        array[j] = len_pr[j];
#else
        array[j] = tbl_pr[j];
#endif
#endif
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
#if DEBUG
#if DEBUG_MATCHES
        array[j] = mch_pr[j];
#elif DEBUG_LENGTH
        array[j] = len_pr[j];
#else
        array[j] = tbl_pr[j];
#endif
#endif
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

        const int * const restrict matrow0 = matrix[s1[i-1+0]];
        const int * const restrict matrow1 = matrix[s1[i-1+1]];
        const int * const restrict matrow2 = matrix[s1[i-1+2]];
        const int * const restrict matrow3 = matrix[s1[i-1+3]];
        const int * const restrict matrow4 = matrix[s1[i-1+4]];
        const int * const restrict matrow5 = matrix[s1[i-1+5]];
        const int * const restrict matrow6 = matrix[s1[i-1+6]];
        const int * const restrict matrow7 = matrix[s1[i-1+7]];

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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength,7);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wscore, 7);
#endif
#endif

        /* j = 1 */
        j = 1;
        /* this block never changes, so let's not copy-and-paste */
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wscore, 6);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wscore, 1);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wscore, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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
#if DEBUG
#if DEBUG_MATCHES
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
#elif DEBUG_LENGTH
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
#else
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
#endif
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

#if DEBUG
    print_array2(array, s1, s1Len, s2, s2Len, 14);
    delete [] array;
#endif
    delete [] s1;
    delete [] s2;
    return DP_t(score, match, length);
}


static int sw_sse8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int score = NEG_INF;
    __m128i vScore = _mm_setzero_si128();
    __m128i vZero = _mm_setzero_si128();
#if DEBUG
    printf("array length (s2Len(=%d)+14+1)*(s1Len(=%d)+7+1) = %d\n",
            s2Len, s1Len, (s2Len+14+1)*(s1Len+7+1));
    int *array = new int[(s2Len+14+1)*(s1Len+7+1)];
    init_vect((s2Len+14+1)*(s1Len+7+1), array, -99);
#endif

    /* dummy padding */
    for (int j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }
    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
#if DEBUG
    array[7] = tbl_pr[7];
#endif
    /* first row */
    for (int j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }
    /* dummy padding */
    for (int j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
#if DEBUG
        array[j] = tbl_pr[j];
#endif
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (int i=1; i<=s1Len; i+=8) {
        int j;
        __m128i NWscore = _mm_set1_epi16(NEG_INF);
        __m128i Nscore  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl    = _mm_set1_epi16(NEG_INF);
        __m128i vDel    = _mm_set1_epi16(NEG_INF);
        __m128i vIns    = _mm_set1_epi16(NEG_INF);
        __m128i vZero   = _mm_set1_epi16(0);
        __m128i vMat;

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, 0, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wscore, 7);
#endif

        /* j = 1 */
        j = 1;
        NWscore = Nscore;
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);
        vDel    = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1],s2[j-1]),
            0,
            0,
            0,
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wscore, 6);
#endif

        /* j = 2 */
        j = 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1+0],s2[j-1-0]),
            BLOSUM(s1[i-1+1],s2[j-1-1]),
            0,
            0,
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif

        /* j = 3 */
        j = 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1+0],s2[j-1-0]),
            BLOSUM(s1[i-1+1],s2[j-1-1]),
            BLOSUM(s1[i-1+2],s2[j-1-2]),
            0,
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif

        /* j = 4 */
        j = 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1+0],s2[j-1-0]),
            BLOSUM(s1[i-1+1],s2[j-1-1]),
            BLOSUM(s1[i-1+2],s2[j-1-2]),
            BLOSUM(s1[i-1+3],s2[j-1-3]),
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif

        /* j = 5 */
        j = 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1+0],s2[j-1-0]),
            BLOSUM(s1[i-1+1],s2[j-1-1]),
            BLOSUM(s1[i-1+2],s2[j-1-2]),
            BLOSUM(s1[i-1+3],s2[j-1-3]),
            BLOSUM(s1[i-1+4],s2[j-1-4]),
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif

        /* j = 6 */
        j = 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1+0],s2[j-1-0]),
            BLOSUM(s1[i-1+1],s2[j-1-1]),
            BLOSUM(s1[i-1+2],s2[j-1-2]),
            BLOSUM(s1[i-1+3],s2[j-1-3]),
            BLOSUM(s1[i-1+4],s2[j-1-4]),
            BLOSUM(s1[i-1+5],s2[j-1-5]),
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wscore, 1);
#endif

        /* j = 7 */
        j = 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM(s1[i-1+0],s2[j-1-0]),
            BLOSUM(s1[i-1+1],s2[j-1-1]),
            BLOSUM(s1[i-1+2],s2[j-1-2]),
            BLOSUM(s1[i-1+3],s2[j-1-3]),
            BLOSUM(s1[i-1+4],s2[j-1-4]),
            BLOSUM(s1[i-1+5],s2[j-1-5]),
            BLOSUM(s1[i-1+6],s2[j-1-6]),
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = 0;
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wscore, 0);
#endif

        for (j=8; j<=s2Len; ++j) {
            NWscore = Nscore;
            Nscore = vshift16(Wscore, tbl_pr[j+7]);
            vDel = vshift16(vDel, del_pr[j+7]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(Nscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(Wscore,vOpen),
                    _mm_sub_epi16(vIns,vGap));
            vMat = _mm_set_epi16(
                    BLOSUM(s1[i-1+0],s2[j-1-0]),
                    BLOSUM(s1[i-1+1],s2[j-1-1]),
                    BLOSUM(s1[i-1+2],s2[j-1-2]),
                    BLOSUM(s1[i-1+3],s2[j-1-3]),
                    BLOSUM(s1[i-1+4],s2[j-1-4]),
                    BLOSUM(s1[i-1+5],s2[j-1-5]),
                    BLOSUM(s1[i-1+6],s2[j-1-6]),
                    BLOSUM(s1[i-1+7],s2[j-1-7])
                    );
            vTbl = _mm_add_epi16(NWscore, vMat);
            vTbl = _mm_max_epi16(vTbl, vDel);
            vTbl = _mm_max_epi16(vTbl, vIns);
            vTbl = _mm_max_epi16(vTbl, vZero);
            vScore = _mm_max_epi16(vScore, vTbl);
            Wscore = vTbl;
            tbl_pr[j] = EXTRACT(vTbl, 0);
            del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
            array[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        }
        if (i+0 == s1Len) {
            int tmp = EXTRACT(vTbl, 7);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 1 */
        j = s2Len + 1;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                BLOSUM(s1[i-1+1],s2[j-1-1]),
                BLOSUM(s1[i-1+2],s2[j-1-2]),
                BLOSUM(s1[i-1+3],s2[j-1-3]),
                BLOSUM(s1[i-1+4],s2[j-1-4]),
                BLOSUM(s1[i-1+5],s2[j-1-5]),
                BLOSUM(s1[i-1+6],s2[j-1-6]),
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+1 == s1Len) {
            int tmp = EXTRACT(vTbl, 6);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 2 */
        j = s2Len + 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                BLOSUM(s1[i-1+2],s2[j-1-2]),
                BLOSUM(s1[i-1+3],s2[j-1-3]),
                BLOSUM(s1[i-1+4],s2[j-1-4]),
                BLOSUM(s1[i-1+5],s2[j-1-5]),
                BLOSUM(s1[i-1+6],s2[j-1-6]),
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+2 == s1Len) {
            int tmp = EXTRACT(vTbl, 5);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 3 */
        j = s2Len + 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                BLOSUM(s1[i-1+3],s2[j-1-3]),
                BLOSUM(s1[i-1+4],s2[j-1-4]),
                BLOSUM(s1[i-1+5],s2[j-1-5]),
                BLOSUM(s1[i-1+6],s2[j-1-6]),
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+3 == s1Len) {
            int tmp = EXTRACT(vTbl, 4);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 4 */
        j = s2Len + 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                BLOSUM(s1[i-1+4],s2[j-1-4]),
                BLOSUM(s1[i-1+5],s2[j-1-5]),
                BLOSUM(s1[i-1+6],s2[j-1-6]),
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+4 == s1Len) {
            int tmp = EXTRACT(vTbl, 3);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 5 */
        j = s2Len + 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                BLOSUM(s1[i-1+5],s2[j-1-5]),
                BLOSUM(s1[i-1+6],s2[j-1-6]),
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+5 == s1Len) {
            int tmp = EXTRACT(vTbl, 2);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 6 */
        j = s2Len + 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM(s1[i-1+6],s2[j-1-6]),
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+6 == s1Len) {
            int tmp = EXTRACT(vTbl, 1);
            score = MAX(score,tmp);
        }

        /* j = s2Len + 7 */
        j = s2Len + 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#if DEBUG
        array[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+7 == s1Len) {
            int tmp = EXTRACT(vTbl, 0);
            score = MAX(score,tmp);
        }
    }

#if DEBUG
    print_array2(array, s1, s1Len, s2, s2Len, 14);
    delete [] array;
#endif
    max8(score, vScore);
    return score;
}


typedef union {
    __m128i m;
    int v[4];
} v4__m128i;


/* shift given vector v and insert val */
static inline __m128i vshift32(const __m128i &v, int val)
{
    v4__m128i ret;
    ret.m = _mm_srli_si128(v, 4);
    ret.v[3] = val;
    return ret.m;
}


/* sse2 doesn't have 32 bit pairwise max */
static inline __m128i vmax32(const __m128i &left, const __m128i &right)
{
    v4__m128i vleft;
    v4__m128i vright;
    v4__m128i ret;
    vleft.m = left;
    vright.m = right;
    ret.v[0] = MAX(vleft.v[0],vright.v[0]);
    ret.v[1] = MAX(vleft.v[1],vright.v[1]);
    ret.v[2] = MAX(vleft.v[2],vright.v[2]);
    ret.v[3] = MAX(vleft.v[3],vright.v[3]);
    return ret.m;
}


/* sse2 doesn't have horizontal max */
static inline int hmax32(const __m128i &m)
{
    int ret;
    v4__m128i tmp;
    tmp.m = m;
    ret = tmp.v[0];
    ret = MAX(ret, tmp.v[1]);
    ret = MAX(ret, tmp.v[2]);
    ret = MAX(ret, tmp.v[3]);
    return ret;
}


static int sw_sse4_pad(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    __m128i Vscore = _mm_setzero_si128();
    __m128i Vzero = _mm_setzero_si128();
    __m128i Vgap_ext = _mm_set1_epi32(gap_ext);
    __m128i Vgap_open = _mm_set1_epi32(gap_open);
    init_vect(lena+6, nogap, 0);
    init_vect(lena+6, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+6)*(lenb(=%d)+3) = %d\n",
            lena, lenb, (lena+6)*(lenb+3));
    int *array = new int[(lena+6)*(lenb+3)];
    init_vect((lena+6)*(lenb+3), array, -9);
#endif
    for (i=0; i<lenb; i+=4) {
        __m128i Vtmp1;
        __m128i Vtmp2;
        __m128i Vmatrix;
        __m128i Va_gap = _mm_sub_epi32(Vzero, Vgap_open);
        v4__m128i Vb_gap;
        Vb_gap.m = _mm_set_epi32(b_gap[2], b_gap[1], b_gap[0], 0);
        v4__m128i Vnogap;
        Vnogap.m = _mm_set_epi32(nogap[2], nogap[1], nogap[0], 0);
        __m128i Vlast_nogap = Vzero;
        __m128i Vprev_nogap = Vzero;
        for (j=0; j<lena+3; ++j) {
            Vb_gap.m = vshift32(Vb_gap.m, b_gap[j+3]);
            Vnogap.m = vshift32(Vnogap.m, nogap[j+3]);
            Vtmp1 = _mm_sub_epi32(Vlast_nogap, Vgap_open);
            Vtmp2 = _mm_sub_epi32(Va_gap, Vgap_ext);
            Va_gap = vmax32(Vtmp1, Vtmp2);
            Vtmp1 = _mm_sub_epi32(Vnogap.m, Vgap_open);
            Vtmp2 = _mm_sub_epi32(Vb_gap.m, Vgap_ext);
            Vb_gap.m = vmax32(Vtmp1, Vtmp2);
            Vmatrix = _mm_set_epi32(
                BLOSUM(seqA[j-0],seqB[i]),
                BLOSUM(seqA[j-1],seqB[i+1]),
                BLOSUM(seqA[j-2],seqB[i+2]),
                BLOSUM(seqA[j-3],seqB[i+3])
            );
            Vtmp1 = _mm_add_epi32(Vprev_nogap, Vmatrix);
            Vlast_nogap = vmax32(Vtmp1, Vzero);
            Vlast_nogap = vmax32(Vlast_nogap, Va_gap);
            Vlast_nogap = vmax32(Vlast_nogap, Vb_gap.m);
            Vprev_nogap = Vnogap.m;
            Vnogap.m = Vlast_nogap;
            b_gap[j] = Vb_gap.v[0];
            nogap[j] = Vnogap.v[0];
#if DEBUG
            array[(i+0)*(lena+6) + j+3] = Vnogap.v[3];
            array[(i+1)*(lena+6) + j+2] = Vnogap.v[2];
            array[(i+2)*(lena+6) + j+1] = Vnogap.v[1];
            array[(i+3)*(lena+6) + j+0] = Vnogap.v[0];
#endif
            Vscore = vmax32(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 6);
    delete [] array;
#endif
    return hmax32(Vscore);
}


static int sw_sse8_pad(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    int score;
    __m128i Vscore = _mm_setzero_si128();
    __m128i Vzero = _mm_setzero_si128();
    __m128i Vgap_ext = _mm_set1_epi16(gap_ext);
    __m128i Vgap_open = _mm_set1_epi16(gap_open);
    init_vect(lena+14, nogap, 0);
    init_vect(lena+14, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+14)*(lenb(=%d)+7) = %d\n",
            lena, lenb, (lena+14)*(lenb+7));
    int *array = new int[(lena+14)*(lenb+7)];
    init_vect((lena+14)*(lenb+7), array, -9);
#endif
    for (i=0; i<lenb; i+=8) {
        __m128i Vtmp1;
        __m128i Vtmp2;
        __m128i Vmatrix;
        __m128i Va_gap = _mm_sub_epi16(Vzero, Vgap_open);
        __m128i Vb_gap = _mm_set_epi16(b_gap[6], b_gap[5], b_gap[4], b_gap[3],
                                       b_gap[2], b_gap[1], b_gap[0], 0);
        __m128i Vnogap = _mm_set_epi16(nogap[6], nogap[5], nogap[4], nogap[3],
                                       nogap[2], nogap[1], nogap[0], 0);
        __m128i Vlast_nogap = Vzero;
        __m128i Vprev_nogap = Vzero;
        for (j=0; j<lena+7; ++j) {
            Vb_gap = vshift16(Vb_gap, b_gap[j+7]);
            Vnogap = vshift16(Vnogap, nogap[j+7]);
            Vtmp1 = _mm_sub_epi16(Vlast_nogap, Vgap_open);
            Vtmp2 = _mm_sub_epi16(Va_gap, Vgap_ext);
            Va_gap = _mm_max_epi16(Vtmp1, Vtmp2);
            Vtmp1 = _mm_sub_epi16(Vnogap, Vgap_open);
            Vtmp2 = _mm_sub_epi16(Vb_gap, Vgap_ext);
            Vb_gap = _mm_max_epi16(Vtmp1, Vtmp2);
            Vmatrix = _mm_set_epi16(
                BLOSUM(seqA[j-0],seqB[i+0]),
                BLOSUM(seqA[j-1],seqB[i+1]),
                BLOSUM(seqA[j-2],seqB[i+2]),
                BLOSUM(seqA[j-3],seqB[i+3]),
                BLOSUM(seqA[j-4],seqB[i+4]),
                BLOSUM(seqA[j-5],seqB[i+5]),
                BLOSUM(seqA[j-6],seqB[i+6]),
                BLOSUM(seqA[j-7],seqB[i+7])
            );
            Vtmp1 = _mm_add_epi16(Vprev_nogap, Vmatrix);
            Vlast_nogap = _mm_max_epi16(Vtmp1, Vzero);
            Vlast_nogap = _mm_max_epi16(Vlast_nogap, Va_gap);
            Vlast_nogap = _mm_max_epi16(Vlast_nogap, Vb_gap);
            Vprev_nogap = Vnogap;
            Vnogap = Vlast_nogap;
            b_gap[j] = EXTRACT(Vb_gap, 0);
            nogap[j] = EXTRACT(Vnogap, 0);
#if DEBUG
            array[(i+0)*(lena+14) + j+7] = EXTRACT(Vnogap, 7);
            array[(i+1)*(lena+14) + j+6] = EXTRACT(Vnogap, 6);
            array[(i+2)*(lena+14) + j+5] = EXTRACT(Vnogap, 5);
            array[(i+3)*(lena+14) + j+4] = EXTRACT(Vnogap, 4);
            array[(i+4)*(lena+14) + j+3] = EXTRACT(Vnogap, 3);
            array[(i+5)*(lena+14) + j+2] = EXTRACT(Vnogap, 2);
            array[(i+6)*(lena+14) + j+1] = EXTRACT(Vnogap, 1);
            array[(i+7)*(lena+14) + j+0] = EXTRACT(Vnogap, 0);
#endif
            Vscore = _mm_max_epi16(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 14);
    delete [] array;
#endif
    max8(score, Vscore);
    return score;
}


int main(int argc, char **argv)
{
#if SHORT_TEST || DEBUG
    const char *seqA = "$$$$$$$$$$$$$$$$SLPSMRADSFTKELMEKISS$$$$$$$$$$$$$$$$";
    const char *seqB = "MTNKICIYAISKNEEKFV$$$$$$$$$$$$$$$$";
#else
    const char *seqA = "$$$$$$$$$$$$$$$$"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "SLPSMRADSFTKELMEKISS"
                       "$$$$$$$$$$$$$$$$";
    const char *seqB = "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "MTNKICIYAISKNEEKFV"
                       "$$$$$$$$$$$$$$$$";
#endif
    const int lena = strlen(seqA) - 32;
    const int lenb = strlen(seqB) - 16;
    const int lenmax = MAX(lena,lenb)+32;
    int *nogap = new int[lenmax];
    int *b_gap = new int[lenmax];
    int *matches = new int[lenmax];
    int *length = new int[lenmax];
    int score;
    DP_t stats;
    unsigned long long timer;
#if DEBUG
    size_t limit = 1;
#else
    size_t limit = 10000;
#endif
    size_t i;
    
    timer_init();
    ::std::cout << timer_name() << " timer" << ::std::endl;

    ::std::cout << "alg\t\ttime\tscore\tmatches\tlength" << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "nw orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "nw_sse8\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = nw_stats(&seqA[16], lena, seqB, lenb, 10, 1, blosum62,
                nogap, b_gap, matches, length);
    }
    timer = timer_end(timer);
    ::std::cout << "nw stat\t"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = nw_stats_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62,
                nogap, b_gap, matches, length);
    }
    timer = timer_end(timer);
    ::std::cout << "nw stat sse"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sg orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sg_sse8\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sg_stats(&seqA[16], lena, seqB, lenb, 10, 1, blosum62,
                nogap, b_gap, matches, length);
    }
    timer = timer_end(timer);
    ::std::cout << "sg stat\t"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sg_stats_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62,
                nogap, b_gap, matches, length);
    }
    timer = timer_end(timer);
    ::std::cout << "sg stat sse"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_woz(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_woz\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_sse4_pad(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_sse4_pad\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_sse8_pad(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_sse8_pad\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_sse8\t\t" << timer/limit << "\t" << score << ::std::endl;
#endif

    delete [] nogap;
    delete [] b_gap;
    delete [] matches;
    delete [] length;

    return 0;
}
