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

#define DEBUG 1
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


static inline void init_vect(int len, int *vec, int val)
{
    int i;
    for (i=0; i<len; ++i) {
        vec[i] = val;
    }
}


#if DEBUG
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
    printf("%s={%3d,%3d,%3d,%3d,%3d,%3d,%3d,%3d}", name,
            int(tmp.v[0]), int(tmp.v[1]), int(tmp.v[2]), int(tmp.v[3]),
            int(tmp.v[4]), int(tmp.v[5]), int(tmp.v[6]), int(tmp.v[7]));
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

    __m128i vTemp;

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
        vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
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
    DP_t max;

    /* array to record the largest score of each reference position */
    uint16_t* maxColumn = (uint16_t*) calloc(s2Len, 2);

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvMStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvMLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvLStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvLLoad = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxScore = vZero;

    /* Trace the highest score till the previous column. */
    __m128i vMaxMark = vZero;

    __m128i vTemp;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        int32_t cmp;
        __m128i e;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = vZero;
        __m128i Nmatch = vZero;
        __m128i Nlength = vZero;
        __m128i Cmatch = vZero;
        __m128i Clength = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* load left and upper left match counts */
        __m128i Wmatch = _mm_load_si128(pvMStore);
        __m128i NWmatch = _mm_slli_si128(pvMStore[segLen - 1], 2);

        /* load left and upper left length counts */
        __m128i Wlength = _mm_load_si128(pvLStore);
        __m128i NWlength = _mm_slli_si128(pvLStore[segLen - 1], 2);

        /* vMaxColumn is used to record the max values of column j. */
        __m128i vMaxColumn = vZero;

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

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

            /* Save vM and vL values */
            __m128i case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,e));
            __m128i case2not = _mm_cmplt_epi16(vF,e);
            __m128i case2 = _mm_andnot_si128(case2not,case1not);
            __m128i case3 = _mm_and_si128(case1not,case2not);
            Cmatch = _mm_andnot_si128(case1not,
                    _mm_add_epi16(NWmatch, _mm_and_si128(
                            // fix me _mm_cmpeq_epi16(vs1,vs2),vOne)));
                            _mm_cmpeq_epi16(vOne,vOne),vOne)));
            Clength= _mm_andnot_si128(case1not,
                    _mm_add_epi16(NWlength, vOne));
            Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));
            Clength= _mm_or_si128(Clength,_mm_and_si128(case2,
                        _mm_add_epi16(Nlength, vOne)));
            Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));
            Clength= _mm_or_si128(Clength,_mm_and_si128(case3,
                        _mm_add_epi16(Wlength, vOne)));
            __m128i mask = _mm_cmplt_epi16(vH,vOne);
            vH = _mm_andnot_si128(mask, vH);
            Cmatch = _mm_andnot_si128(mask, Cmatch);
            Clength= _mm_andnot_si128(mask, Clength);
            mask = _mm_cmpgt_epi16(vH,vScore);
            vScore = _mm_andnot_si128(mask, vScore);
            vScore = _mm_or_si128(vScore, _mm_and_si128(mask, vH));
            vMatch = _mm_andnot_si128(mask, vMatch);
            vMatch = _mm_or_si128(vMatch, _mm_and_si128(mask, Cmatch));
            vLength= _mm_andnot_si128(mask, vLength);
            vLength= _mm_or_si128(vLength,_mm_and_si128(mask, Clength));

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
        vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
        cmp = _mm_movemask_epi8(vTemp);
        if (cmp != 0xffff) {
            uint16_t temp;
            vMaxMark = vMaxScore;
            max8(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (temp > max.score) {
                max.score = temp;
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


int main(int /*argc*/, char ** /*argv*/)
{
#if SHORT_TEST || DEBUG
    const char *seqA = "SLPSMRADSFTKELMEKISS";
    const char *seqB = "MTNKICIYAISKNEEKFV";
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

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "nw orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sg orig\t\t" << timer/limit << "\t" << score << ::std::endl;
#endif

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw stat\t"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;

    return 0;
}
