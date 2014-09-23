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

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
					(m) = _mm_extract_epi16((vm), 0)

static const int16_t NEG_INF = SHRT_MIN / 2;
//static const int16_t NEG_INF = -50;

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
                    //*t++ = j>=s1Len ? 0 :
                    //    matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    //*s++ = j>=s1Len ? 0 :
                    //    nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = (nt == MAP_BLOSUM_[(unsigned char)s1[j]]);
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
                //*h = j>=s1Len ? NEG_INF : -open-gap*(segNum*segLen+i);
                //*e = j>=s1Len ? NEG_INF : *h-open;
                *h = -open-gap*(segNum*segLen+i);
                *e = *h-open;
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
            //7*segLen < s1Len ? -open-open-7*segLen*gap : NEG_INF,
            //6*segLen < s1Len ? -open-open-6*segLen*gap : NEG_INF,
            //5*segLen < s1Len ? -open-open-5*segLen*gap : NEG_INF,
            //4*segLen < s1Len ? -open-open-4*segLen*gap : NEG_INF,
            //3*segLen < s1Len ? -open-open-3*segLen*gap : NEG_INF,
            //2*segLen < s1Len ? -open-open-2*segLen*gap : NEG_INF,
            //1*segLen < s1Len ? -open-open-1*segLen*gap : NEG_INF,
            //0*segLen < s1Len ? -open-open-0*segLen*gap : NEG_INF);
            -open-open-7*segLen*gap,
            -open-open-6*segLen*gap,
            -open-open-5*segLen*gap,
            -open-open-4*segLen*gap,
            -open-open-3*segLen*gap,
            -open-open-2*segLen*gap,
            -open-open-1*segLen*gap,
            -open-open-0*segLen*gap);
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
            vF = _mm_insert_epi16(vF, boundary[j+1]-open, 0);
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
                    //*t++ = j>=s1Len ? 0 :
                    //    matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
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
                //*h = j>=s1Len ? NEG_INF : -open-gap*(segNum*segLen+i);
                //*e = j>=s1Len ? NEG_INF : *h-open;
                *h = -open-gap*(segNum*segLen+i);
                *e = *h-open;
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
            //7*segLen < s1Len ? -open-open-7*segLen*gap : NEG_INF,
            //6*segLen < s1Len ? -open-open-6*segLen*gap : NEG_INF,
            //5*segLen < s1Len ? -open-open-5*segLen*gap : NEG_INF,
            //4*segLen < s1Len ? -open-open-4*segLen*gap : NEG_INF,
            //3*segLen < s1Len ? -open-open-3*segLen*gap : NEG_INF,
            //2*segLen < s1Len ? -open-open-2*segLen*gap : NEG_INF,
            //1*segLen < s1Len ? -open-open-1*segLen*gap : NEG_INF,
            //0*segLen < s1Len ? -open-open-0*segLen*gap : NEG_INF);
            -open-open-7*segLen*gap,
            -open-open-6*segLen*gap,
            -open-open-5*segLen*gap,
            -open-open-4*segLen*gap,
            -open-open-3*segLen*gap,
            -open-open-2*segLen*gap,
            -open-open-1*segLen*gap,
            -open-open-0*segLen*gap);
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
            vF = _mm_insert_epi16(vF, boundary[j+1]-open, 0);
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
                    //*t++ = j>=s1Len ? 0 :
                    //    matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    //*s++ = j>=s1Len ? 0 :
                    //    nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = (nt == MAP_BLOSUM_[(unsigned char)s1[j]]);
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

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int max = NEG_INF;

    /* Define 16 byte 0 vector. */
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        for (int32_t nt=0; nt<n; ++nt) {
            for (int32_t i=0; i<segLen; ++i) {
                int32_t j = i;
                for (int32_t segNum=0; segNum<8; ++segNum) {
                    //*t++ = j>=s1Len ? 0 :
                    //    matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize E */
    {
        int16_t *e = (int16_t*)pvE;
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
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = initialF;

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
            vF = _mm_insert_epi16(vF, -open, 0);
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
            /* extract last value from the column */
            vH = _mm_load_si128(pvHStore + offset);
            for (int32_t k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 2);
            }
            /* max of last value in each column */
            int16_t tmp = (signed short) _mm_extract_epi16 (vH, 7);
            if (tmp > max) {
                max = tmp;
            }
        }
    }

    /* max of last column */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vMaxLastColH = vNegInf;
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
        __m128i cond = _mm_cmplt_epi16(vH, vMaxLastColH);
        vMaxLastColH = _mm_or_si128(_mm_andnot_si128(cond, vH),
                                    _mm_and_si128(cond, vMaxLastColH));
        vQIndex = _mm_adds_epi16(vQIndex, vOne);
    }

    /* max in vec */
    for (int32_t j=0; j<8; ++j) {
        int16_t value = (signed short) _mm_extract_epi16(vMaxLastColH, j);
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

#if DEBUG
    printf("array length (s1Len(=%d)+0)*s2Len(=%d) = %d\n",
            s1Len, s2Len, (s1Len+0)*s2Len);
    int *array = new int[segLen*8*s2Len];
    init_vect(segLen*8, array, -9);
#endif

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
                    //*t++ = j>=s1Len ? 0 :
                    //    matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    //*s++ = j>=s1Len ? 0 :
                    //    nt == MAP_BLOSUM_[(unsigned char)s1[j]] ? 1 : 0;
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = (nt == MAP_BLOSUM_[(unsigned char)s1[j]]);
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
        pv = pvELoad;
        pvELoad = pvEStore;
        pvEStore = pv;

        __m128i vQIndex = vQIndex_reset;

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
            vH = _mm_max_epi16(vH, vZero);
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

    /* max in vec */
    for (int32_t j=0; j<8; ++j) {
        int16_t value = (signed short) _mm_extract_epi16(vMaxH, j);
        if (value > max.score) {
            max.score = value;
            max.matches = (signed short) _mm_extract_epi16(vMaxM, j);
            max.length = (signed short) _mm_extract_epi16(vMaxL, j);
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

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        for (int32_t nt=0; nt<n; ++nt) {
            for (int32_t i=0; i<segLen; ++i) {
                int32_t j = i;
                for (int32_t segNum=0; segNum<8; ++segNum) {
                    //*t++ = j>=s1Len ? 0 :
                    //    matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize E */
    {
        int16_t *e = (int16_t*)pvE;
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

    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxH = vZero;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        __m128i vE;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
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
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            vH = _mm_max_epi16(vH, vZero);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
#if DEBUG
            arr_store_si128(array, vH, i, segLen, j, s2Len);
#endif

            /* update max vector seen so far */
            __m128i cond_max = _mm_cmpgt_epi16(vH,vMaxH);
            __m128i cond_lmt = _mm_cmplt_epi16(vQIndex,vQLimit);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxH = _mm_andnot_si128(cond_all, vMaxH);
            vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
            vQIndex = _mm_add_epi16(vQIndex, vOne);

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

    /* max in vec */
    for (int32_t j=0; j<8; ++j) {
        int16_t value = (signed short) _mm_extract_epi16(vMaxH, j);
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

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sg striped\t" << timer/limit << "\t" << score << ::std::endl;
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

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw striped\t" << timer/limit << "\t" << score << ::std::endl;
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        stats = sw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw striped stat"
        << "\t" << timer/limit
        << "\t" << stats.score
        << "\t" << stats.matches
        << "\t" << stats.length
        << ::std::endl;
#endif

    return 0;
}
