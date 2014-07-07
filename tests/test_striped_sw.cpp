/* All DP algorithms here implement Farrar's approach of vectorizing along the
 * column/query sequence. */
#include "config.h"

#include <emmintrin.h>

#include <stdint.h>

#include <climits>
#include <cstring>
#include <iostream>

#include "blosum/blosum62.h"
#include "timer.h"

using namespace ::std;

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


int sw(
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
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[s1[j]]];
                j += segLen;
            }
        }
    }

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

        /* Right part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[s2[j]] * segLen;

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

    return max;
}


int sg(
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
                *t++ = j>=s1Len ? 0 : matrix[nt*n + MAP_BLOSUM_[s1[j]]];
                j += segLen;
            }
        }
    }

    /* the max alignment score */
    uint16_t max = SHRT_MIN/4;

    /* array to record the largest score of each reference position */
    uint16_t* maxColumn = (uint16_t*) calloc(s2Len, 2);

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();

    /* Define 16 byte SHRT_MIN vector. */
    __m128i vShortMin = _mm_set1_epi16(SHRT_MIN/4);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* initialize E */
    t = (int16_t*)pvE;
    for (int32_t i=0; i<segLen; ++i) {
        int32_t j = i;
        for (int32_t segNum=0; segNum<8; ++segNum) {
            *t++ = j>=s1Len ? 0 : -open - j*gap;
            j += segLen;
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxScore = vShortMin;

    /* Trace the highest score till the previous column. */
    __m128i vMaxMark = vShortMin;

    __m128i vTemp;

    /* outer loop over database sequence */
    for (int32_t j=0; j<s2Len; ++j) {
        int32_t cmp;
        __m128i e;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = vShortMin;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* vMaxColumn is used to record the max values of column j. */
        __m128i vMaxColumn = vShortMin;

        /* Right part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[s2[j]] * segLen;

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

    return max;
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
        score = sw(&seqA[16], lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sw orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(&seqA[16], lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    ::std::cout << "sg orig\t\t" << timer/limit << "\t" << score << ::std::endl;

    return 0;
}
