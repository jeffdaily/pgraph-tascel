/* All DP algorithms here implement Farrar's approach of vectorizing along the
 * column/query sequence. */
#include "config.h"

#include <emmintrin.h>

#include <stdint.h>

#include <climits>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "blosum/blosum62.h"
#include "timer.h"

using namespace ::std;

#define DEBUG 0
#define DEBUG_MATCHES 0
#define DEBUG_LENGTH 0
#define SHORT_TEST 0

#define MAX(a,b) ((a)>(b)?(a):(b))

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
					(m) = _mm_extract_epi16((vm), 0)

//static const int16_t NEG_INF = SHRT_MIN / 2;
static const int16_t NEG_INF = -50;

struct DP_t {
    int score;
    int matches;
    int length;
    DP_t() : score(0), matches(0), length(0) {}
    DP_t(int score, int matches, int length)
        : score(score), matches(matches), length(length) {}
};


#if DEBUG
static inline void init_vect(int len, int *vec, int val)
{
    int i;
    for (i=0; i<len; ++i) {
        vec[i] = val;
    }
}


static inline int extract(const __m128i &m, const int &pos)
{
    union {
        __m128i m;
        int16_t v[8];
    } tmp;
    tmp.m = m;
    return tmp.v[pos];
}


static void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = extract(vH, 0);
    array[(1*seglen+t)*dlen + d] = extract(vH, 1);
    array[(2*seglen+t)*dlen + d] = extract(vH, 2);
    array[(3*seglen+t)*dlen + d] = extract(vH, 3);
    array[(4*seglen+t)*dlen + d] = extract(vH, 4);
    array[(5*seglen+t)*dlen + d] = extract(vH, 5);
    array[(6*seglen+t)*dlen + d] = extract(vH, 6);
    array[(7*seglen+t)*dlen + d] = extract(vH, 7);
}


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
#if 0
            int(tmp.v[0]), int(tmp.v[1]), int(tmp.v[2]), int(tmp.v[3]),
            int(tmp.v[4]), int(tmp.v[5]), int(tmp.v[6]), int(tmp.v[7])
#else
            extract(m,0), extract(m,1), extract(m,2), extract(m,3),
            extract(m,4), extract(m,5), extract(m,6), extract(m,7)
#endif
            );
}
#endif


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


static int sw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = (int16_t*)vProfile;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    /* the max alignment score */
    uint16_t max = 0;

    /* array to record the largest score of each reference position */
    uint16_t* maxColumn = (uint16_t*) calloc(s2Len, 2);

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxScore = vZero;

    /* Trace the highest score till the previous column. */
    __m128i vMaxMark = vZero;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        int32_t cmp;
        __m128i e;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* vMaxColumn is used to record the max values of column j. */
        __m128i vMaxColumn = vZero;

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));

            /* Get max from vH, vE and vF. */
            e = _mm_load_si128(pvE + i);
            vH = _mm_max_epi16(vH, e);
            vH = _mm_max_epi16(vH, vF);
            vMaxColumn = _mm_max_epi16(vMaxColumn, vH);

            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
#if DEBUG
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            /* saturation arithmetic, result >= 0 */
            vH = _mm_subs_epu16(vH, vGapO);
            e = _mm_subs_epu16(e, vGapE);
            e = _mm_max_epi16(e, vH);
            _mm_store_si128(pvE + i, e);

            /* Update vF value. */
            vF = _mm_subs_epu16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128 (vF, 2);
            for (int32_t i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH, vF);
                _mm_store_si128(pvHStore + i, vH);
#if DEBUG
                arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epu16(vH, vGapO);
                vF = _mm_subs_epu16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
            }
        }

end:
        vMaxScore = _mm_max_epi16(vMaxScore, vMaxColumn);
        __m128i vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
        cmp = _mm_movemask_epi8(vTemp);
        if (cmp != 0xffff) {
            uint16_t temp;
            vMaxMark = vMaxScore;
            max8(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (temp > max) {
                max = temp;
                for (int32_t i=0; i<segLen; ++i) {
                    pvHmax[i] = pvHStore[i];
                }
            }
        }

        /* Record the max score of current column. */
        max8(maxColumn[j], vMaxColumn);
    }

    free(vProfile);
    free(pvHmax);
    free(pvE);
    free(pvHLoad);
    free(pvHStore);

    free(maxColumn);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return max;
}


static DP_t sw_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* vProfileS = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = NULL;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    t = (int16_t*)vProfile;
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }
    t = (int16_t*)vProfileS;
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 :
                    nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    /* the max alignment score */
    DP_t max;

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEM = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEL = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* for max calculation we don't want to include padded cells */
    __m128i vQLimit = _mm_set1_epi16(s1Len);
    __m128i vQIndex_reset = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxH = vZero;
    __m128i vMaxM = vZero;
    __m128i vMaxL = vZero;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vEM;
        __m128i vEL;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = vZero;
        __m128i vFM = vZero;
        __m128i vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        __m128i vHM = _mm_slli_si128(pvHMStore[segLen - 1], 2);
        __m128i vHL = _mm_slli_si128(pvHLStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        const __m128i* vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
        pv = pvHLLoad;
        pvHLLoad = pvHLStore;
        pvHLStore = pv;

        __m128i vQIndex = vQIndex_reset;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            __m128i case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            __m128i case2not = _mm_cmplt_epi16(vF,vE);
            __m128i case2 = _mm_andnot_si128(case2not,case1not);
            __m128i case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            vH = _mm_max_epi16(vH, vZero);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            __m128i mask = _mm_cmplt_epi16(vH,vOne);
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_add_epi16(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            vHM = _mm_andnot_si128(mask, vHM);
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_add_epi16(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_add_epi16(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_add_epi16(vEL, vOne)));
            vHL = _mm_andnot_si128(mask, vHL);
            _mm_store_si128(pvHLStore + i, vHL);

#if DEBUG
#if DEBUG_MATCHES
            arr_store_si128(array, vHM, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
            arr_store_si128(array, vHL, i, segLen, j, s2Len);
#else
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif

            /* update max vector seen so far */
            __m128i cond_max = _mm_cmpgt_epi16(vH,vMaxH);
            __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxH = _mm_andnot_si128(cond_all, vMaxH);
            vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
            vMaxM = _mm_andnot_si128(cond_all, vMaxM);
            vMaxM = _mm_or_si128(vMaxM, _mm_and_si128(cond_all, vHM));
            vMaxL = _mm_andnot_si128(cond_all, vMaxL);
            vMaxL = _mm_or_si128(vMaxL, _mm_and_si128(cond_all, vHL));
            vQIndex = _mm_add_epi16(vQIndex, vOne);

            /* Update vE value. */
            /* saturation arithmetic, result >= 0 */
            vH = _mm_subs_epu16(vH, vGapO);
            vE = _mm_subs_epu16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvE + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epu16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128 (vF, 2);
            vFM = _mm_slli_si128 (vFM, 2);
            vFL = _mm_slli_si128 (vFL, 2);
            vQIndex = vQIndex_reset;
            for (int32_t i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vHM = _mm_load_si128(pvHMStore + i);
                vHL = _mm_load_si128(pvHLStore + i);
                __m128i cond = _mm_cmpgt_epi16(vF,vH);
                vH = _mm_andnot_si128(cond, vH);
                vH = _mm_or_si128(vH, _mm_and_si128(cond, vF));
                _mm_store_si128(pvHStore + i, vH);
                vHM = _mm_andnot_si128(cond, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(cond, vFM));
                _mm_store_si128(pvHMStore + i, vHM);
                vHL = _mm_andnot_si128(cond, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(cond, vFL));
                _mm_store_si128(pvHLStore + i, vHL);
#if DEBUG
#if DEBUG_MATCHES
                arr_store_si128(array, vHM, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
                arr_store_si128(array, vHL, i, segLen, j, s2Len);
#else
                arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif
                /* update max vector seen so far */
                __m128i cond_max = _mm_cmpgt_epi16(vH,vMaxH);
                __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxH = _mm_andnot_si128(cond_all, vMaxH);
                vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
                vMaxM = _mm_andnot_si128(cond_all, vMaxM);
                vMaxM = _mm_or_si128(vMaxM, _mm_and_si128(cond_all, vHM));
                vMaxL = _mm_andnot_si128(cond_all, vMaxL);
                vMaxL = _mm_or_si128(vMaxL, _mm_and_si128(cond_all, vHL));
                vQIndex = _mm_add_epi16(vQIndex, vOne);

                vH = _mm_subs_epu16(vH, vGapO);
                vF = _mm_subs_epu16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vFM = vHM;
                vFL = vHL;
            }
        }
end:
        {
        }
    }

    /* max in vec */
    for (int32_t j=0; j<8; ++j) {
        int value = _mm_extract_epi16(vMaxH, j);
        if (value > max.score) {
            max.score = value;
            max.matches = _mm_extract_epi16(vMaxM, j);
            max.length = _mm_extract_epi16(vMaxL, j);
        }
    }

    free(vProfile);
    free(vProfileS);
    free(pvHStore);
    free(pvHLoad);
    free(pvHMStore);
    free(pvHMLoad);
    free(pvHLStore);
    free(pvHLLoad);
    free(pvE);
    free(pvEM);
    free(pvEL);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return max;
}


static int sw_mine_ORIG(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = (int16_t*)vProfile;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    /* the max alignment score */
    uint16_t max = 0;

    /* array to record the largest score of each reference position */
    uint16_t* maxColumn = (uint16_t*) calloc(s2Len, 2);

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxScore = vZero;

    /* Trace the highest score till the previous column. */
    __m128i vMaxMark = vZero;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        int32_t cmp;

        /* Initialize F value to 0.
         * Any errors to vH values will be corrected in the Lazy_F loop. */
        __m128i vF = vZero;
        __m128i vH_i1 = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        /* aka H_(i-1,j-1) */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* vMaxColumn is used to record the max values of column j. */
        __m128i vMaxColumn = vZero;

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            /* calculate vE */
            __m128i vH_j1 = _mm_load_si128(pvHLoad + i);
            __m128i vE_j1 = _mm_load_si128(pvE + i);
            vH_j1 = _mm_subs_epi16(vH_j1, vGapO);
            vE_j1 = _mm_subs_epi16(vE_j1, vGapE);
            __m128i vE = _mm_max_epi16(vE_j1, vH_j1);
            _mm_store_si128(pvE + i, vE);

            /* calculate vF */
            vF = _mm_subs_epi16(vF, vGapE);
            vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
            vF = _mm_max_epi16(vF, vH_i1);

            /* calculate new vH */
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            vH = _mm_max_epi16(vH, vZero);
            vMaxColumn = _mm_max_epi16(vMaxColumn, vH);
            _mm_store_si128(pvHStore + i, vH);
#if DEBUG
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
            /* keep latest vH for next inner loop */
            vH_i1 = vH;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 2);
            vH_i1 = _mm_slli_si128(vH_i1, 2);
            for (int32_t i=0; i<segLen; ++i) {
                vF = _mm_subs_epi16(vF,vGapE);
                vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
                vF = _mm_max_epi16(vF,vH_i1);
                vH = _mm_load_si128(pvHStore + i);
                __m128i cond = _mm_cmpgt_epi16(vF,vH);
                if (!_mm_movemask_epi8(cond)) {
                    goto end;
                }
                vH_i1 = _mm_max_epi16(vF,vH);
                _mm_store_si128(pvHStore + i, vH_i1);
#if DEBUG
                arr_store_si128(array, vH_i1, i, segLen, j, s2Len);
#endif
            }
        }

end:
        vMaxScore = _mm_max_epi16(vMaxScore, vMaxColumn);
        __m128i vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
        cmp = _mm_movemask_epi8(vTemp);
        if (cmp != 0xffff) {
            uint16_t temp;
            vMaxMark = vMaxScore;
            max8(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (temp > max) {
                max = temp;
                for (int32_t i=0; i<segLen; ++i) {
                    pvHmax[i] = pvHStore[i];
                }
            }
        }

        /* Record the max score of current column. */
        max8(maxColumn[j], vMaxColumn);
    }

    free(vProfile);
    free(pvHmax);
    free(pvE);
    free(pvHLoad);
    free(pvHStore);

    free(maxColumn);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return max;
}


static int sw_mine(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = NULL;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    t = (int16_t*)vProfile;
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    /* the max alignment score */
    int max = NEG_INF;

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* for max calculation we don't want to include padded cells */
    __m128i vQLimit = _mm_set1_epi16(s1Len);
    __m128i vQIndex_reset = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);

    /* Define 16 byte NEG_INF vector. */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);

    /* vMaxColumn is used to record the max values of column j. */
    __m128i vMaxColumn = vNegInf;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        /* Initialize F value to 0.
         * Any errors to vH values will be corrected in the Lazy_F loop. */
        __m128i vF = vZero;
        __m128i vH_i1 = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        /* aka H_(i-1,j-1) */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        __m128i vQIndex = vQIndex_reset;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            /* calculate vE */
            __m128i vH_j1 = _mm_load_si128(pvHLoad + i);
            __m128i vE_j1 = _mm_load_si128(pvE + i);
            vH_j1 = _mm_subs_epi16(vH_j1, vGapO);
            vE_j1 = _mm_subs_epi16(vE_j1, vGapE);
            __m128i vE = _mm_max_epi16(vE_j1, vH_j1);
            _mm_store_si128(pvE + i, vE);

            /* calculate vF */
            vF = _mm_subs_epi16(vF, vGapE);
            vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
            vF = _mm_max_epi16(vF, vH_i1);

            /* calculate new vH */
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            vH = _mm_max_epi16(vH, vZero);
            _mm_store_si128(pvHStore + i, vH);

#if DEBUG
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            __m128i cond_max = _mm_cmpgt_epi16(vH,vMaxColumn);
            __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxColumn = _mm_andnot_si128(cond_all, vMaxColumn);
            vMaxColumn = _mm_or_si128(vMaxColumn,
                    _mm_and_si128(cond_all, vH));
            vQIndex = _mm_add_epi16(vQIndex, vOne);

            /* keep latest vH for next inner loop */
            vH_i1 = vH;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 2);
            vH_i1 = _mm_slli_si128(vH_i1, 2);
            vQIndex = vQIndex_reset;
            for (int32_t i=0; i<segLen; ++i) {
                vF = _mm_subs_epi16(vF,vGapE);
                vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
                vF = _mm_max_epi16(vF,vH_i1);
                vH = _mm_load_si128(pvHStore + i);
                __m128i cond = _mm_cmpgt_epi16(vF,vH);
                if (!_mm_movemask_epi8(cond)) {
                    goto end;
                }
                vH_i1 = _mm_max_epi16(vF,vH);
                _mm_store_si128(pvHStore + i, vH_i1);
#if DEBUG
                arr_store_si128(array, vH_i1, i, segLen, j, s2Len);
#endif
                /* update max vector seen so far */
                __m128i cond_max = _mm_cmpgt_epi16(vH_i1,vMaxColumn);
                __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxColumn = _mm_andnot_si128(cond_all, vMaxColumn);
                vMaxColumn = _mm_or_si128(vMaxColumn,
                        _mm_and_si128(cond_all, vH_i1));
                vQIndex = _mm_add_epi16(vQIndex, vOne);
            }
        }
end:
        {
        }
    }

    for (int32_t j=0; j<8; ++j) {
        int value = _mm_extract_epi16(vMaxColumn, j);
        if (value > max) {
            max = value;
        }
    }

    free(vProfile);
    free(pvHStore);
    free(pvHLoad);
    free(pvE);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return max;
}


static DP_t sw_mine_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* vProfileM= (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = NULL;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    t = (int16_t*)vProfile;
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }
    /* Generate query profile rearrange query sequence & calculate the 
     * basic match/mismatch as a 1/0 */
    t = (int16_t*)vProfileM;
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 :
                    nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    /* the max alignment score */
    DP_t max;

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvMStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvMLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvLStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvLLoad = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* for max calculation we don't want to include padded cells */
    __m128i vQLimit = _mm_set1_epi16(s1Len);
    __m128i vQIndex_reset = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);

    /* Define 16 byte NEG_INF vector. */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);

    /* vMaxColumn is used to record the max values of column j. */
    __m128i vMaxColumn = vNegInf;
    __m128i vMaxColumnM= vNegInf;
    __m128i vMaxColumnL= vNegInf;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        /* Initialize F value to 0.
         * Any errors to vH values will be corrected in the Lazy_F loop. */
        __m128i vF = vZero;
        __m128i vH_i1 = vZero;
        __m128i Nmatch = vZero;
        __m128i Nlength = vZero;
        __m128i Cmatch = vZero;
        __m128i Clength = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        /* aka H_(i-1,j-1) */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        __m128i NWmatch = _mm_slli_si128(pvMStore[segLen - 1], 2);
        __m128i NWlength = _mm_slli_si128(pvLStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        const __m128i* vPM= vProfileM+ MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* Swap the 2 M buffers. */
        pv = pvMLoad;
        pvMLoad = pvMStore;
        pvMStore = pv;

        /* Swap the 2 L buffers. */
        pv = pvLLoad;
        pvLLoad = pvLStore;
        pvLStore = pv;

        __m128i vQIndex = vQIndex_reset;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            /* calculate vE */
            __m128i vH_j1 = _mm_load_si128(pvHLoad + i);
            __m128i vE_j1 = _mm_load_si128(pvE + i);
            vH_j1 = _mm_subs_epi16(vH_j1, vGapO);
            vE_j1 = _mm_subs_epi16(vE_j1, vGapE);
            __m128i vE = _mm_max_epi16(vE_j1, vH_j1);
            _mm_store_si128(pvE + i, vE);

            /* calculate vF */
            vF = _mm_subs_epi16(vF, vGapE);
            vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
            vF = _mm_max_epi16(vF, vH_i1);

            /* calculate new vH */
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            __m128i case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            __m128i case2not = _mm_cmplt_epi16(vF,vE);
            __m128i case2 = _mm_andnot_si128(case2not,case1not);
            __m128i case3 = _mm_and_si128(case1not,case2not);

            /* finish vH calc */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            vH = _mm_max_epi16(vH, vZero);
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM and vL */
            __m128i Wmatch = _mm_load_si128(pvMLoad + i);
            __m128i Wlength = _mm_load_si128(pvLLoad + i);
            Cmatch = _mm_andnot_si128(case1not,
                    _mm_add_epi16(NWmatch, _mm_load_si128(vPM + i)));
            Clength= _mm_andnot_si128(case1not,
                    _mm_add_epi16(NWlength, vOne));
            Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));
            Clength= _mm_or_si128(Clength,_mm_and_si128(case2,
                        _mm_add_epi16(Nlength, vOne)));
            Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));
            Clength= _mm_or_si128(Clength,_mm_and_si128(case3,
                        _mm_add_epi16(Wlength, vOne)));
            __m128i mask = _mm_cmplt_epi16(vH,vOne);
            Cmatch = _mm_andnot_si128(mask, Cmatch);
            Clength= _mm_andnot_si128(mask, Clength);
            _mm_store_si128(pvMStore + i, Cmatch);
            _mm_store_si128(pvLStore + i, Clength);
#if DEBUG
#if DEBUG_MATCHES
            arr_store_si128(array, Cmatch, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
            arr_store_si128(array, Clength, i, segLen, j, s2Len);
#else
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif
            /* update max vector seen so far */
            __m128i cond_max = _mm_cmpgt_epi16(vH,vMaxColumn);
            __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
            __m128i cond = _mm_and_si128(cond_max, cond_lmt);
            vMaxColumn = _mm_andnot_si128(cond, vMaxColumn);
            vMaxColumn = _mm_or_si128(vMaxColumn,
                    _mm_and_si128(cond, vH));
            vMaxColumnM= _mm_andnot_si128(cond, vMaxColumnM);
            vMaxColumnM= _mm_or_si128(vMaxColumnM,
                    _mm_and_si128(cond, Cmatch));
            vMaxColumnL= _mm_andnot_si128(cond, vMaxColumnL);
            vMaxColumnL= _mm_or_si128(vMaxColumnL,
                    _mm_and_si128(cond, Clength));
            vQIndex = _mm_add_epi16(vQIndex, vOne);

            /* keep latest vH for next inner loop */
            vH_i1 = vH;
            Nmatch = Cmatch;
            Nlength= Clength;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            NWmatch = Wmatch;
            NWlength = Wlength;
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 2);
            vH_i1 = _mm_slli_si128(vH_i1, 2);
            Nmatch = _mm_slli_si128(Nmatch, 2);
            Nlength= _mm_slli_si128(Nlength, 2);
            vQIndex = vQIndex_reset;
            for (int32_t i=0; i<segLen; ++i) {
                vF = _mm_subs_epi16(vF,vGapE);
                vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
                vF = _mm_max_epi16(vF,vH_i1);
                vH = _mm_load_si128(pvHStore + i);
                __m128i cond = _mm_cmpgt_epi16(vF,vH);
                if (!_mm_movemask_epi8(cond)) {
                    goto end;
                }
                vH_i1 = _mm_max_epi16(vF,vH);
                _mm_store_si128(pvHStore + i, vH_i1);

                Cmatch = _mm_load_si128(pvMStore + i);
                Cmatch = _mm_andnot_si128(cond, Cmatch);
                Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(cond,
                        _mm_add_epi16(Nmatch, _mm_load_si128(vPM + i))));
                _mm_store_si128(pvMStore + i, Cmatch);
                Nmatch = Cmatch;

                Clength= _mm_load_si128(pvLStore + i);
                Clength= _mm_andnot_si128(cond, Clength);
                Clength= _mm_or_si128(Clength, _mm_and_si128(cond,
                        _mm_add_epi16(Nlength, vOne)));
                _mm_store_si128(pvLStore + i, Clength);
                Nlength= Clength;
#if DEBUG
#if DEBUG_MATCHES
                arr_store_si128(array, Cmatch, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
                arr_store_si128(array, Clength, i, segLen, j, s2Len);
#else
                arr_store_si128(array, vH_i1, i, segLen, j, s2Len);
#endif
#endif
                /* update max vector seen so far */
                __m128i cond_max = _mm_cmpgt_epi16(vH_i1,vMaxColumn);
                __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxColumn = _mm_andnot_si128(cond_all, vMaxColumn);
                vMaxColumn = _mm_or_si128(vMaxColumn,
                        _mm_and_si128(cond_all, vH_i1));
                vMaxColumnM= _mm_andnot_si128(cond_all, vMaxColumnM);
                vMaxColumnM= _mm_or_si128(vMaxColumnM,
                        _mm_and_si128(cond_all, Cmatch));
                vMaxColumnL= _mm_andnot_si128(cond_all, vMaxColumnL);
                vMaxColumnL= _mm_or_si128(vMaxColumnL,
                        _mm_and_si128(cond_all, Clength));
                vQIndex = _mm_add_epi16(vQIndex, vOne);
            }
        }
end:
        {
        }
    }

    for (int32_t j=0; j<8; ++j) {
        int value = _mm_extract_epi16(vMaxColumn, j);
        if (value > max.score) {
            max.score = value;
            max.matches =_mm_extract_epi16(vMaxColumnM, j);
            max.length = _mm_extract_epi16(vMaxColumnL, j);
        }
    }

    free(vProfile);
    free(vProfileM);
    free(pvHStore);
    free(pvHLoad);
    free(pvE);
    free(pvMStore);
    free(pvMLoad);
    free(pvLStore);
    free(pvLLoad);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return max;
}


static int sg(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = (int16_t*)vProfile;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int16_t maxLastRow = NEG_INF;

    /* Define 16 byte NEG_INF vector. */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* initialize E */
    t = (int16_t*)pvE;
    for (int32_t i=0; i<segLen; ++i) {
        for (int32_t segNum=0; segNum<8; ++segNum) {
            *t++ = NEG_INF;
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i e;
        /* Initialize F value to -INF.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegInf;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));

            /* Get max from vH, vE and vF. */
            e = _mm_load_si128(pvE + i);
            vH = _mm_max_epi16(vH, e);
            vH = _mm_max_epi16(vH, vF);

            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            e = _mm_subs_epi16(e, vGapE);
            e = _mm_max_epi16(e, vH);
            _mm_store_si128(pvE + i, e);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128 (vF, 2);
            /* first vF value must be -INF */
            vF = _mm_insert_epi16(vF, NEG_INF, 0);
            for (int32_t i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH, vF);
                _mm_store_si128(pvHStore + i, vH);
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
            }
        }
end:
        /* extract last value from the column */
        vH = _mm_load_si128(pvHStore + offset);
        for (int32_t k=0; k<position; ++k) {
            vH = _mm_slli_si128 (vH, 2);
        }
        /* max of last value in each column */
        int16_t tmp = (signed short) _mm_extract_epi16 (vH, 7);
        maxLastRow = MAX(maxLastRow, tmp);
#if DEBUG
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvHStore + i);
            arr_store_si128(array, vH, i, segLen, j, s2Len);
        }
#endif
    }

    /* max of last column */
    __m128i vMaxLastCol = vNegInf;
    __m128i vQIndex = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);
    __m128i vQLimit = _mm_set1_epi16(s1Len);
    __m128i vOne = _mm_set1_epi16(1);

    for (int32_t i=0; i<segLen; ++i) {
        /* load the last stored values */
        __m128i vH = _mm_load_si128(pvHStore + i);
        /* mask off the values that were padded */
        __m128i vMask = _mm_cmplt_epi16(vQIndex, vQLimit);
        vH = _mm_or_si128(_mm_and_si128(vMask, vH),
                          _mm_andnot_si128(vMask, vNegInf));
        vMaxLastCol = _mm_max_epi16(vMaxLastCol, vH);
        vQIndex = _mm_adds_epi16(vQIndex, vOne);
    }
    int16_t maxLastCol;
    max8(maxLastCol, vMaxLastCol);

    free(vProfile);
    free(pvE);
    free(pvHLoad);
    free(pvHStore);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return MAX(maxLastCol,maxLastRow);
}


static int sg_mine(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = NULL;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    t = (int16_t*)vProfile;
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int16_t maxLastRow = NEG_INF;

    /* Define 16 byte 0 vector. */
    //__m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    /* Define 16 byte NEG_INF vector. */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* initialize E */
    t = (int16_t*)pvE;
    for (int32_t i=0; i<segLen; ++i) {
        for (int32_t segNum=0; segNum<8; ++segNum) {
            *t++ = NEG_INF;
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        /* Initialize F value to 0.
         * Any errors to vH values will be corrected in the Lazy_F loop. */
        __m128i vF = vNegInf;
        __m128i vH_i1 = vNegInf;

        /* load final segment of pvHStore and shift left by 2 bytes */
        /* aka H_(i-1,j-1) */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            /* calculate vE */
            __m128i vH_j1 = _mm_load_si128(pvHLoad + i);
            __m128i vE_j1 = _mm_load_si128(pvE + i);
            vH_j1 = _mm_subs_epi16(vH_j1, vGapO);
            vE_j1 = _mm_subs_epi16(vE_j1, vGapE);
            __m128i vE = _mm_max_epi16(vE_j1, vH_j1);
            _mm_store_si128(pvE + i, vE);

            /* calculate vF */
            vF = _mm_subs_epi16(vF, vGapE);
            vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
            vF = _mm_max_epi16(vF, vH_i1);

            /* calculate new vH */
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            _mm_store_si128(pvHStore + i, vH);

#if DEBUG
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
            /* keep latest vH for next inner loop */
            vH_i1 = vH;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 2);
            /* first vF value must be -INF */
            vF = _mm_insert_epi16(vF, NEG_INF, 0);
            vH_i1 = _mm_slli_si128(vH_i1, 2);
            for (int32_t i=0; i<segLen; ++i) {
                vF = _mm_subs_epi16(vF,vGapE);
                vH_i1 = _mm_subs_epi16(vH_i1, vGapO);
                vF = _mm_max_epi16(vF,vH_i1);
                vH = _mm_load_si128(pvHStore + i);
                __m128i cond = _mm_cmpgt_epi16(vF,vH);
                if (!_mm_movemask_epi8(cond)) {
                    goto end;
                }
                vH_i1 = _mm_max_epi16(vF,vH);
                _mm_store_si128(pvHStore + i, vH_i1);
#if DEBUG
                arr_store_si128(array, vH_i1, i, segLen, j, s2Len);
#endif
            }
        }
end:
        {
            /* extract last value from the column */
            vH = _mm_load_si128(pvHStore + offset);
            for (int32_t k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 2);
            }
            /* max of last value in each column */
            int16_t tmp = (signed short) _mm_extract_epi16 (vH, 7);
            maxLastRow = MAX(maxLastRow, tmp);
#if DEBUG
            for (int32_t i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                arr_store_si128(array, vH, i, segLen, j, s2Len);
            }
#endif
        }
    }

    /* max of last column */
    __m128i vMaxLastCol = vNegInf;
    __m128i vQIndex = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);
    __m128i vQLimit = _mm_set1_epi16(s1Len);

    for (int32_t i=0; i<segLen; ++i) {
        /* load the last stored values */
        __m128i vH = _mm_load_si128(pvHStore + i);
        /* mask off the values that were padded */
        __m128i vMask = _mm_cmplt_epi16(vQIndex, vQLimit);
        vH = _mm_or_si128(_mm_and_si128(vMask, vH),
                          _mm_andnot_si128(vMask, vNegInf));
        vMaxLastCol = _mm_max_epi16(vMaxLastCol, vH);
        vQIndex = _mm_adds_epi16(vQIndex, vOne);
    }
    int16_t maxLastCol;
    max8(maxLastCol, vMaxLastCol);

    free(vProfile);
    free(pvHStore);
    free(pvHLoad);
    free(pvE);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return MAX(maxLastCol,maxLastRow);
}


#if 0
static int nw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    int16_t* t = (int16_t*)vProfile;

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    for (int32_t nt=0; nt<n; ++nt) {
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                j += segLen;
            }
        }
    }

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int16_t last_value = NEG_INF;

    /* Define 16 byte NEG_INF vector. */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    int16_t* boundary = (int16_t*) calloc(s2Len, sizeof(int16_t));

    /* initialize H */
    {
        t = (int16_t*)pvHStore;
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *t++ = j>=s1Len ? NEG_INF : -open-gap*(segNum*segLen+i);
                j += segLen;
            }
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (int32_t i=1; i<s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
        }
    }

    /* initialize E */
    t = (int16_t*)pvE;
    for (int32_t i=0; i<segLen; ++i) {
        for (int32_t segNum=0; segNum<8; ++segNum) {
            *t++ = NEG_INF;
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i e;
        /* Initialize F value to -INF.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegInf;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* insert upper boundary condition */
        vH = _mm_insert_epi16(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));

            /* Get max from vH, vE and vF. */
            e = _mm_load_si128(pvE + i);
            vH = _mm_max_epi16(vH, e);
            vH = _mm_max_epi16(vH, vF);

            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            e = _mm_subs_epi16(e, vGapE);
            e = _mm_max_epi16(e, vH);
            _mm_store_si128(pvE + i, e);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128 (vF, 2);
            /* first vF value must be -INF */
            vF = _mm_insert_epi16(vF, NEG_INF, 0);
            for (int32_t i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH, vF);
                _mm_store_si128(pvHStore + i, vH);
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
            }
        }
end:
        {
        }
#if DEBUG
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvHStore + i);
            arr_store_si128(array, vH, i, segLen, j, s2Len);
        }
#endif
    }

    /* extract last value from the last column */
    __m128i vH = _mm_load_si128(pvHStore + offset);
    for (int32_t k=0; k<position; ++k) {
        vH = _mm_slli_si128 (vH, 2);
    }
    last_value = (signed short) _mm_extract_epi16 (vH, 7);

    free(vProfile);
    free(pvE);
    free(pvHLoad);
    free(pvHStore);
    free(boundary);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return last_value;
}
#endif


static DP_t nw_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* vProfileS = (__m128i*)malloc(n * segLen * sizeof(__m128i));

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    DP_t last_value;

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvELoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEM = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEL = (__m128i*) calloc(segLen, sizeof(__m128i));
    int16_t* boundary = (int16_t*) calloc(s2Len+1, sizeof(int16_t));

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        int16_t *s = (int16_t*)vProfileS;
        for (int32_t nt=0; nt<n; ++nt) {
            for (int32_t i=0; i<segLen; ++i) {
                int32_t j = i;
                for (int32_t segNum=0; segNum<8; ++segNum) {
                    *t++ = j>=s1Len ? 0 :
                        matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = j>=s1Len ? 0 :
                        nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                    j += segLen;
                }
            }
        }
    }

    /* initialize H and E */
    {
        int16_t *h = (int16_t*)pvHStore;
        int16_t *e = (int16_t*)pvEStore;
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *h = j>=s1Len ? NEG_INF : -open-gap*(segNum*segLen+i);
                *e = j>=s1Len ? NEG_INF : MAX(*h-open,NEG_INF-gap);
                j += segLen;
                ++h;
                ++e;
            }
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (int32_t i=1; i<=s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    __m128i initialF = _mm_set_epi16(
            7*segLen < s1Len ? -open-open-7*segLen*gap : NEG_INF,
            6*segLen < s1Len ? -open-open-6*segLen*gap : NEG_INF,
            5*segLen < s1Len ? -open-open-5*segLen*gap : NEG_INF,
            4*segLen < s1Len ? -open-open-4*segLen*gap : NEG_INF,
            3*segLen < s1Len ? -open-open-3*segLen*gap : NEG_INF,
            2*segLen < s1Len ? -open-open-2*segLen*gap : NEG_INF,
            1*segLen < s1Len ? -open-open-1*segLen*gap : NEG_INF,
            0*segLen < s1Len ? -open-open-0*segLen*gap : NEG_INF);
    initialF = _mm_adds_epi16(initialF, vGapE);

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vEM;
        __m128i vEL;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        initialF = _mm_subs_epi16(initialF, vGapE);
        __m128i vF = initialF;
        __m128i vFM = vZero;
        __m128i vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        __m128i vHM = _mm_slli_si128(pvHMStore[segLen - 1], 2);
        __m128i vHL = _mm_slli_si128(pvHLStore[segLen - 1], 2);

        /* insert upper boundary condition */
        vH = _mm_insert_epi16(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        const __m128i* vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
        pv = pvHLLoad;
        pvHLLoad = pvHLStore;
        pvHLStore = pv;
        pv = pvELoad;
        pvELoad = pvEStore;
        pvEStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            __m128i case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            __m128i case2not = _mm_cmplt_epi16(vF,vE);
            __m128i case2 = _mm_andnot_si128(case2not,case1not);
            __m128i case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_add_epi16(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_add_epi16(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_add_epi16(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_add_epi16(vEL, vOne)));
            _mm_store_si128(pvHLStore + i, vHL);

#if DEBUG
#if DEBUG_MATCHES
            arr_store_si128(array, vHM, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
            arr_store_si128(array, vHL, i, segLen, j, s2Len);
#else
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            vE = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 2);
            vHp = _mm_insert_epi16(vHp, boundary[j], 0);
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, MAX(boundary[j+1]-open,NEG_INF-gap), 0);
            vFM = _mm_slli_si128(vFM, 2);
            vFL = _mm_slli_si128(vFL, 2);
            for (int32_t i=0; i<segLen; ++i) {
                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_adds_epi16(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                __m128i case1not = _mm_or_si128(
                        _mm_cmplt_epi16(vHp,vF),_mm_cmplt_epi16(vHp,vE));
                __m128i case2not = _mm_cmplt_epi16(vF,vE);
                __m128i case2 = _mm_andnot_si128(case2not,case1not);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_andnot_si128(case2, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_andnot_si128(case2, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                            _mm_add_epi16(vFL,vOne)));
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#if DEBUG
#if DEBUG_MATCHES
                arr_store_si128(array, vHM, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
                arr_store_si128(array, vHL, i, segLen, j, s2Len);
#else
                arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm_max_epi16(vF, vH);
                vFM = vHM;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    __m128i vH = _mm_load_si128(pvHStore + offset);
    __m128i vHM = _mm_load_si128(pvHMStore + offset);
    __m128i vHL = _mm_load_si128(pvHLStore + offset);
    for (int32_t k=0; k<position; ++k) {
        vH = _mm_slli_si128 (vH, 2);
        vHM = _mm_slli_si128 (vHM, 2);
        vHL = _mm_slli_si128 (vHL, 2);
    }
    last_value.score = (signed short) _mm_extract_epi16 (vH, 7);
    last_value.matches = (signed short) _mm_extract_epi16 (vHM, 7);
    last_value.length = (signed short) _mm_extract_epi16 (vHL, 7);

    free(vProfile);
    free(vProfileS);
    free(pvHStore);
    free(pvHLoad);
    free(pvHMStore);
    free(pvHMLoad);
    free(pvHLStore);
    free(pvHLLoad);
    free(pvEStore);
    free(pvELoad);
    free(pvEM);
    free(pvEL);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return last_value;
}


static int nw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int last_value;

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    int16_t* boundary = (int16_t*) calloc(s2Len+1, sizeof(int16_t));

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        for (int32_t nt=0; nt<n; ++nt) {
            for (int32_t i=0; i<segLen; ++i) {
                int32_t j = i;
                for (int32_t segNum=0; segNum<8; ++segNum) {
                    *t++ = j>=s1Len ? 0 :
                        matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize H and E */
    {
        int16_t *h = (int16_t*)pvHStore;
        int16_t *e = (int16_t*)pvE;
        for (int32_t i=0; i<segLen; ++i) {
            int32_t j = i;
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *h = j>=s1Len ? NEG_INF : -open-gap*(segNum*segLen+i);
                *e = j>=s1Len ? NEG_INF : MAX(*h-open,NEG_INF-gap);
                j += segLen;
                ++h;
                ++e;
            }
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (int32_t i=1; i<=s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    __m128i initialF = _mm_set_epi16(
            7*segLen < s1Len ? -open-open-7*segLen*gap : NEG_INF,
            6*segLen < s1Len ? -open-open-6*segLen*gap : NEG_INF,
            5*segLen < s1Len ? -open-open-5*segLen*gap : NEG_INF,
            4*segLen < s1Len ? -open-open-4*segLen*gap : NEG_INF,
            3*segLen < s1Len ? -open-open-3*segLen*gap : NEG_INF,
            2*segLen < s1Len ? -open-open-2*segLen*gap : NEG_INF,
            1*segLen < s1Len ? -open-open-1*segLen*gap : NEG_INF,
            0*segLen < s1Len ? -open-open-0*segLen*gap : NEG_INF);
    initialF = _mm_adds_epi16(initialF, vGapE);

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i vE;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        initialF = _mm_subs_epi16(initialF, vGapE);
        __m128i vF = initialF;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* insert upper boundary condition */
        vH = _mm_insert_epi16(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
#if DEBUG
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            vE = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, MAX(boundary[j+1]-open,NEG_INF-gap), 0);
            for (int32_t i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#if DEBUG
                arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm_max_epi16(vF, vH);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    __m128i vH = _mm_load_si128(pvHStore + offset);
    for (int32_t k=0; k<position; ++k) {
        vH = _mm_slli_si128 (vH, 2);
    }
    last_value = (signed short) _mm_extract_epi16 (vH, 7);

    free(vProfile);
    free(pvHStore);
    free(pvHLoad);
    free(pvE);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return last_value;
}


static DP_t sg_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix)
{
    /* compute 'query' profile; s1 is the query */
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* vProfileS = (__m128i*)malloc(n * segLen * sizeof(__m128i));

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    DP_t max;
    max.score = NEG_INF;

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvELoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEM = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEL = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        int16_t *s = (int16_t*)vProfileS;
        for (int32_t nt=0; nt<n; ++nt) {
            for (int32_t i=0; i<segLen; ++i) {
                int32_t j = i;
                for (int32_t segNum=0; segNum<8; ++segNum) {
                    *t++ = j>=s1Len ? 0 :
                        matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = j>=s1Len ? 0 :
                        nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                    j += segLen;
                }
            }
        }
    }

    /* initialize E */
    {
        int16_t *e = (int16_t*)pvEStore;
        for (int32_t i=0; i<segLen; ++i) {
            for (int32_t segNum=0; segNum<8; ++segNum) {
                *e++ = -open;
            }
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    __m128i initialF = _mm_set_epi16(
            -open-7*segLen*gap,
            -open-6*segLen*gap,
            -open-5*segLen*gap,
            -open-4*segLen*gap,
            -open-3*segLen*gap,
            -open-2*segLen*gap,
            -open-1*segLen*gap,
            -open-0*segLen*gap);

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vEM;
        __m128i vEL;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = initialF;
        __m128i vFM = vZero;
        __m128i vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        __m128i vHM = _mm_slli_si128(pvHMStore[segLen - 1], 2);
        __m128i vHL = _mm_slli_si128(pvHLStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        const __m128i* vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
        pv = pvHLLoad;
        pvHLLoad = pvHLStore;
        pvHLStore = pv;
        pv = pvELoad;
        pvELoad = pvEStore;
        pvEStore = pv;

        /* inner loop to process the query sequence */
        for (int32_t i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            __m128i case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            __m128i case2not = _mm_cmplt_epi16(vF,vE);
            __m128i case2 = _mm_andnot_si128(case2not,case1not);
            __m128i case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_add_epi16(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_add_epi16(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_add_epi16(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_add_epi16(vEL, vOne)));
            _mm_store_si128(pvHLStore + i, vHL);

#if DEBUG
#if DEBUG_MATCHES
            arr_store_si128(array, vHM, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
            arr_store_si128(array, vHL, i, segLen, j, s2Len);
#else
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            vE = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (int32_t k=0; k<8; ++k) {
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 2);
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, -open, 0);
            vFM = _mm_slli_si128(vFM, 2);
            vFL = _mm_slli_si128(vFL, 2);
            for (int32_t i=0; i<segLen; ++i) {
                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_adds_epi16(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                __m128i case1not = _mm_or_si128(
                        _mm_cmplt_epi16(vHp,vF),_mm_cmplt_epi16(vHp,vE));
                __m128i case2not = _mm_cmplt_epi16(vF,vE);
                __m128i case2 = _mm_andnot_si128(case2not,case1not);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_andnot_si128(case2, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_andnot_si128(case2, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                            _mm_add_epi16(vFL,vOne)));
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#if DEBUG
#if DEBUG_MATCHES
                arr_store_si128(array, vHM, i, segLen, j, s2Len);
#elif DEBUG_LENGTH
                arr_store_si128(array, vHL, i, segLen, j, s2Len);
#else
                arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif
#endif
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm_max_epi16(vF, vH);
                vFM = vHM;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
            /* extract last value from the column */
            vH = _mm_load_si128(pvHStore + offset);
            vHM = _mm_load_si128(pvHMStore + offset);
            vHL = _mm_load_si128(pvHLStore + offset);
            for (int32_t k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 2);
                vHM = _mm_slli_si128 (vHM, 2);
                vHL = _mm_slli_si128 (vHL, 2);
            }
            /* max of last value in each column */
            int16_t tmp = (signed short) _mm_extract_epi16 (vH, 7);
            if (tmp > max.score) {
                max.score = tmp;
                max.matches = _mm_extract_epi16(vHM, 7);
                max.length = _mm_extract_epi16(vHL, 7);
            }
        }
    }

    /* max of last column */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vMaxLastColH = vNegInf;
    __m128i vMaxLastColHM = vNegInf;
    __m128i vMaxLastColHL = vNegInf;
    __m128i vQIndex = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);
    __m128i vQLimit = _mm_set1_epi16(s1Len);

    for (int32_t i=0; i<segLen; ++i) {
        /* load the last stored values */
        __m128i vH = _mm_load_si128(pvHStore + i);
        __m128i vHM = _mm_load_si128(pvHMStore + i);
        __m128i vHL = _mm_load_si128(pvHLStore + i);
        /* mask off the values that were padded */
        __m128i vMask = _mm_cmplt_epi16(vQIndex, vQLimit);
        vH = _mm_or_si128(_mm_and_si128(vMask, vH),
                          _mm_andnot_si128(vMask, vNegInf));
        vHM = _mm_or_si128(_mm_and_si128(vMask, vHM),
                           _mm_andnot_si128(vMask, vNegInf));
        vHL = _mm_or_si128(_mm_and_si128(vMask, vHL),
                           _mm_andnot_si128(vMask, vNegInf));
        __m128i cond = _mm_cmplt_epi16(vH, vMaxLastColH);
        vMaxLastColH = _mm_or_si128(_mm_andnot_si128(cond, vH),
                                    _mm_and_si128(cond, vMaxLastColH));
        vMaxLastColHM = _mm_or_si128(_mm_andnot_si128(cond, vHM),
                                     _mm_and_si128(cond, vMaxLastColHM));
        vMaxLastColHL = _mm_or_si128(_mm_andnot_si128(cond, vHL),
                                     _mm_and_si128(cond, vMaxLastColHL));
        vQIndex = _mm_adds_epi16(vQIndex, vOne);
    }

    /* max in vec */
    for (int32_t j=0; j<8; ++j) {
        int16_t value = (signed short) _mm_extract_epi16(vMaxLastColH, j);
        if (value > max.score) {
            max.score = value;
            max.matches = _mm_extract_epi16(vMaxLastColHM, j);
            max.length = _mm_extract_epi16(vMaxLastColHL, j);
        }
    }

    free(vProfile);
    free(vProfileS);
    free(pvHStore);
    free(pvHLoad);
    free(pvHMStore);
    free(pvHMLoad);
    free(pvHLStore);
    free(pvHLLoad);
    free(pvEStore);
    free(pvELoad);
    free(pvEM);
    free(pvEL);

#if DEBUG
    print_array(array, s2, s2Len, s1, s1Len, 0);
    delete [] array;
#endif

    return max;
}


int main(int /*argc*/, char ** /*argv*/)
{
#if SHORT_TEST || DEBUG
    //const char *seqA = "SLPSMRADSFTKELMEKISS";
    //const char *seqB = "MTNKICIYAISKNEEKFV";
    const char *seqA = "WTHECK";
    const char *seqB = "BETTERWORK";
#else
    const char *seqA =
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS";
    const char *seqB =
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV";
#endif
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
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

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "nw striped\t" << timer/limit << "\t" << score << ::std::endl;
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = nw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "nw striped stat"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sg orig\t\t" << timer/limit << "\t" << score << ::std::endl;
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sg_stats(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sg striped stat"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_mine(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sg mine\t\t" << timer/limit << "\t" << score << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw orig\t\t" << timer/limit << "\t" << score << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw stat"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_mine(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw mine\t\t" << timer/limit << "\t" << score << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sw_mine_stats(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw mine stat"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;
#endif

    return 0;
}
