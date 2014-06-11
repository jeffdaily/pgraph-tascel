#include "config.h"

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

#define DEBUG 1
#define SHORT_TEST 1

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

#define MAX(a,b) ((a)>(b)?(a):(b))

#define max16(m, vm) (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 8)); \
                      (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 4)); \
                      (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 2)); \
                      (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 1)); \
                      (m) = _mm_extract_epi16((vm), 0)

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
                    (m) = _mm_extract_epi16((vm), 0)

//static const int16_t NEG_INF = SHRT_MIN / 2;
static const int16_t NEG_INF = -50;

static void init_vect(int len, int *vec, int val)
{
    int i;
    for (i=0; i<len; ++i) {
        vec[i] = val;
    }
}


struct Vint2 {
    int v[2];

    Vint2() {
        v[0] = 0;
        v[1] = 0;
    }

    Vint2(const int &val) {
        v[0] = val;
        v[1] = val;
    }

    Vint2(const int &val0, const int &val1) {
        v[0] = val0;
        v[1] = val1;
    }

    Vint2 operator=(const Vint2 &that) {
        v[0] = that.v[0];
        v[1] = that.v[1];
        return *this;
    }

    Vint2 operator-() {
        Vint2 result;
        result.v[0] = -v[0];
        result.v[1] = -v[1];
        return result;
    }

    int max() {
        int result = v[0];
        if (v[1] > result) result = v[1];
        return result;
    }

    friend Vint2 operator+(const Vint2 &left, const Vint2 &right);
    friend Vint2 operator-(const Vint2 &left, const Vint2 &right);
    friend Vint2 vshift(const Vint2 &that, const int &val);
    friend Vint2 vmax(const Vint2 &left, const Vint2 &right);
};

Vint2 operator+(const Vint2 &left, const Vint2 &right) {
    Vint2 result;
    result.v[0] = left.v[0] + right.v[0];
    result.v[1] = left.v[1] + right.v[1];
    return result;
}

Vint2 operator-(const Vint2 &left, const Vint2 &right) {
    Vint2 result;
    result.v[0] = left.v[0] - right.v[0];
    result.v[1] = left.v[1] - right.v[1];
    return result;
}

Vint2 vshift(const Vint2 &that, const int &val) {
    Vint2 result;
    result.v[1] = that.v[0];
    result.v[0] = val;
    return result;
}

Vint2 vmax(const Vint2 &left, const Vint2 &right) {
    Vint2 result;
    result.v[0] = MAX(left.v[0],right.v[0]);
    result.v[1] = MAX(left.v[1],right.v[1]);
    return result;
}


struct Vint4 {
    int v[4];

    Vint4() {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
        v[3] = 0;
    }

    Vint4(const int &val) {
        v[0] = val;
        v[1] = val;
        v[2] = val;
        v[3] = val;
    }

    Vint4(const int &val0, const int &val1, const int &val2, const int &val3) {
        v[0] = val0;
        v[1] = val1;
        v[2] = val2;
        v[3] = val3;
    }

    Vint4 operator=(const Vint4 &that) {
        v[0] = that.v[0];
        v[1] = that.v[1];
        v[2] = that.v[2];
        v[3] = that.v[3];
        return *this;
    }

    Vint4 operator-() {
        Vint4 result;
        result.v[0] = -v[0];
        result.v[1] = -v[1];
        result.v[2] = -v[2];
        result.v[3] = -v[3];
        return result;
    }

    int max() {
        int result = v[0];
        if (v[1] > result) result = v[1];
        if (v[2] > result) result = v[2];
        if (v[3] > result) result = v[3];
        return result;
    }

    friend Vint4 operator+(const Vint4 &left, const Vint4 &right);
    friend Vint4 operator-(const Vint4 &left, const Vint4 &right);
    friend Vint4 vshift(const Vint4 &that, const int &val);
    friend Vint4 vmax(const Vint4 &left, const Vint4 &right);
};

Vint4 operator+(const Vint4 &left, const Vint4 &right) {
    Vint4 result;
    result.v[0] = left.v[0] + right.v[0];
    result.v[1] = left.v[1] + right.v[1];
    result.v[2] = left.v[2] + right.v[2];
    result.v[3] = left.v[3] + right.v[3];
    return result;
}

Vint4 operator-(const Vint4 &left, const Vint4 &right) {
    Vint4 result;
    result.v[0] = left.v[0] - right.v[0];
    result.v[1] = left.v[1] - right.v[1];
    result.v[2] = left.v[2] - right.v[2];
    result.v[3] = left.v[3] - right.v[3];
    return result;
}

Vint4 vshift(const Vint4 &that, const int &val) {
    Vint4 result;
    result.v[3] = that.v[2];
    result.v[2] = that.v[1];
    result.v[1] = that.v[0];
    result.v[0] = val;
    return result;
}

Vint4 vmax(const Vint4 &left, const Vint4 &right) {
    Vint4 result;
    result.v[0] = MAX(left.v[0],right.v[0]);
    result.v[1] = MAX(left.v[1],right.v[1]);
    result.v[2] = MAX(left.v[2],right.v[2]);
    result.v[3] = MAX(left.v[3],right.v[3]);
    return result;
}


struct Vint8 {
    int v[8];

    Vint8() {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
        v[3] = 0;
        v[4] = 0;
        v[5] = 0;
        v[6] = 0;
        v[7] = 0;
    }

    Vint8(const int &val) {
        v[0] = val;
        v[1] = val;
        v[2] = val;
        v[3] = val;
        v[4] = val;
        v[5] = val;
        v[6] = val;
        v[7] = val;
    }

    Vint8(const int &val0, const int &val1, const int &val2, const int &val3,
          const int &val4, const int &val5, const int &val6, const int &val7) {
        v[0] = val0;
        v[1] = val1;
        v[2] = val2;
        v[3] = val3;
        v[4] = val4;
        v[5] = val5;
        v[6] = val6;
        v[7] = val7;
    }

    Vint8 operator=(const Vint8 &that) {
        v[0] = that.v[0];
        v[1] = that.v[1];
        v[2] = that.v[2];
        v[3] = that.v[3];
        v[4] = that.v[4];
        v[5] = that.v[5];
        v[6] = that.v[6];
        v[7] = that.v[7];
        return *this;
    }

    Vint8 operator-() {
        Vint8 result;
        result.v[0] = -v[0];
        result.v[1] = -v[1];
        result.v[2] = -v[2];
        result.v[3] = -v[3];
        result.v[4] = -v[4];
        result.v[5] = -v[5];
        result.v[6] = -v[6];
        result.v[7] = -v[7];
        return result;
    }

    int max() {
        int result = v[0];
        if (v[1] > result) result = v[1];
        if (v[2] > result) result = v[2];
        if (v[3] > result) result = v[3];
        if (v[4] > result) result = v[4];
        if (v[5] > result) result = v[5];
        if (v[6] > result) result = v[6];
        if (v[7] > result) result = v[7];
        return result;
    }

    friend Vint8 operator+(const Vint8 &left, const Vint8 &right);
    friend Vint8 operator-(const Vint8 &left, const Vint8 &right);
    friend Vint8 vshift(const Vint8 &that, const int &val);
    friend Vint8 vmax(const Vint8 &left, const Vint8 &right);
};

Vint8 operator+(const Vint8 &left, const Vint8 &right) {
    Vint8 result;
    result.v[0] = left.v[0] + right.v[0];
    result.v[1] = left.v[1] + right.v[1];
    result.v[2] = left.v[2] + right.v[2];
    result.v[3] = left.v[3] + right.v[3];
    result.v[4] = left.v[4] + right.v[4];
    result.v[5] = left.v[5] + right.v[5];
    result.v[6] = left.v[6] + right.v[6];
    result.v[7] = left.v[7] + right.v[7];
    return result;
}

Vint8 operator-(const Vint8 &left, const Vint8 &right) {
    Vint8 result;
    result.v[0] = left.v[0] - right.v[0];
    result.v[1] = left.v[1] - right.v[1];
    result.v[2] = left.v[2] - right.v[2];
    result.v[3] = left.v[3] - right.v[3];
    result.v[4] = left.v[4] - right.v[4];
    result.v[5] = left.v[5] - right.v[5];
    result.v[6] = left.v[6] - right.v[6];
    result.v[7] = left.v[7] - right.v[7];
    return result;
}

Vint8 vshift(const Vint8 &that, const int &val) {
    Vint8 result;
    result.v[7] = that.v[6];
    result.v[6] = that.v[5];
    result.v[5] = that.v[4];
    result.v[4] = that.v[3];
    result.v[3] = that.v[2];
    result.v[2] = that.v[1];
    result.v[1] = that.v[0];
    result.v[0] = val;
    return result;
}

Vint8 vmax(const Vint8 &left, const Vint8 &right) {
    Vint8 result;
    result.v[0] = MAX(left.v[0],right.v[0]);
    result.v[1] = MAX(left.v[1],right.v[1]);
    result.v[2] = MAX(left.v[2],right.v[2]);
    result.v[3] = MAX(left.v[3],right.v[3]);
    result.v[4] = MAX(left.v[4],right.v[4]);
    result.v[5] = MAX(left.v[5],right.v[5]);
    result.v[6] = MAX(left.v[6],right.v[6]);
    result.v[7] = MAX(left.v[7],right.v[7]);
    return result;
}


struct Vint16 {
    int v[16];

    Vint16() {
        v[ 0] = 0;
        v[ 1] = 0;
        v[ 2] = 0;
        v[ 3] = 0;
        v[ 4] = 0;
        v[ 5] = 0;
        v[ 6] = 0;
        v[ 7] = 0;
        v[ 8] = 0;
        v[ 9] = 0;
        v[10] = 0;
        v[11] = 0;
        v[12] = 0;
        v[13] = 0;
        v[14] = 0;
        v[15] = 0;
    }

    Vint16(const int &val) {
        v[ 0] = val;
        v[ 1] = val;
        v[ 2] = val;
        v[ 3] = val;
        v[ 4] = val;
        v[ 5] = val;
        v[ 6] = val;
        v[ 7] = val;
        v[ 8] = val;
        v[ 9] = val;
        v[10] = val;
        v[11] = val;
        v[12] = val;
        v[13] = val;
        v[14] = val;
        v[15] = val;
    }

    Vint16(const int &val0, const int &val1, const int &val2, const int &val3,
           const int &val4, const int &val5, const int &val6, const int &val7,
           const int &val8, const int &val9, const int &val10, const int &val11,
           const int &val12, const int &val13, const int &val14, const int &val15) {
        v[ 0] = val0;
        v[ 1] = val1;
        v[ 2] = val2;
        v[ 3] = val3;
        v[ 4] = val4;
        v[ 5] = val5;
        v[ 6] = val6;
        v[ 7] = val7;
        v[ 8] = val8;
        v[ 9] = val9;
        v[10] = val10;
        v[11] = val11;
        v[12] = val12;
        v[13] = val13;
        v[14] = val14;
        v[15] = val15;
    }

    Vint16 operator=(const Vint16 &that) {
        v[ 0] = that.v[ 0];
        v[ 1] = that.v[ 1];
        v[ 2] = that.v[ 2];
        v[ 3] = that.v[ 3];
        v[ 4] = that.v[ 4];
        v[ 5] = that.v[ 5];
        v[ 6] = that.v[ 6];
        v[ 7] = that.v[ 7];
        v[ 8] = that.v[ 8];
        v[ 9] = that.v[ 9];
        v[10] = that.v[10];
        v[11] = that.v[11];
        v[12] = that.v[12];
        v[13] = that.v[13];
        v[14] = that.v[14];
        v[15] = that.v[15];
        return *this;
    }

    Vint16 operator-() {
        Vint16 result;
        result.v[ 0] = -v[ 0];
        result.v[ 1] = -v[ 1];
        result.v[ 2] = -v[ 2];
        result.v[ 3] = -v[ 3];
        result.v[ 4] = -v[ 4];
        result.v[ 5] = -v[ 5];
        result.v[ 6] = -v[ 6];
        result.v[ 7] = -v[ 7];
        result.v[ 8] = -v[ 8];
        result.v[ 9] = -v[ 9];
        result.v[10] = -v[10];
        result.v[11] = -v[11];
        result.v[12] = -v[12];
        result.v[13] = -v[13];
        result.v[14] = -v[14];
        result.v[15] = -v[15];
        return result;
    }

    int max() {
        int result = v[0];
        if (v[ 1] > result) result = v[ 1];
        if (v[ 2] > result) result = v[ 2];
        if (v[ 3] > result) result = v[ 3];
        if (v[ 4] > result) result = v[ 4];
        if (v[ 5] > result) result = v[ 5];
        if (v[ 6] > result) result = v[ 6];
        if (v[ 7] > result) result = v[ 7];
        if (v[ 8] > result) result = v[ 8];
        if (v[ 9] > result) result = v[ 9];
        if (v[10] > result) result = v[10];
        if (v[11] > result) result = v[11];
        if (v[12] > result) result = v[12];
        if (v[13] > result) result = v[13];
        if (v[14] > result) result = v[14];
        if (v[15] > result) result = v[15];
        return result;
    }

    friend Vint16 operator+(const Vint16 &left, const Vint16 &right);
    friend Vint16 operator-(const Vint16 &left, const Vint16 &right);
    friend Vint16 vshift(const Vint16 &that, const int &val);
    friend Vint16 vmax(const Vint16 &left, const Vint16 &right);
};

Vint16 operator+(const Vint16 &left, const Vint16 &right) {
    Vint16 result;
    result.v[ 0] = left.v[ 0] + right.v[ 0];
    result.v[ 1] = left.v[ 1] + right.v[ 1];
    result.v[ 2] = left.v[ 2] + right.v[ 2];
    result.v[ 3] = left.v[ 3] + right.v[ 3];
    result.v[ 4] = left.v[ 4] + right.v[ 4];
    result.v[ 5] = left.v[ 5] + right.v[ 5];
    result.v[ 6] = left.v[ 6] + right.v[ 6];
    result.v[ 7] = left.v[ 7] + right.v[ 7];
    result.v[ 8] = left.v[ 8] + right.v[ 8];
    result.v[ 9] = left.v[ 9] + right.v[ 9];
    result.v[10] = left.v[10] + right.v[10];
    result.v[11] = left.v[11] + right.v[11];
    result.v[12] = left.v[12] + right.v[12];
    result.v[13] = left.v[13] + right.v[13];
    result.v[14] = left.v[14] + right.v[14];
    result.v[15] = left.v[15] + right.v[15];
    return result;
}

Vint16 operator-(const Vint16 &left, const Vint16 &right) {
    Vint16 result;
    result.v[ 0] = left.v[ 0] - right.v[ 0];
    result.v[ 1] = left.v[ 1] - right.v[ 1];
    result.v[ 2] = left.v[ 2] - right.v[ 2];
    result.v[ 3] = left.v[ 3] - right.v[ 3];
    result.v[ 4] = left.v[ 4] - right.v[ 4];
    result.v[ 5] = left.v[ 5] - right.v[ 5];
    result.v[ 6] = left.v[ 6] - right.v[ 6];
    result.v[ 7] = left.v[ 7] - right.v[ 7];
    result.v[ 8] = left.v[ 8] - right.v[ 8];
    result.v[ 9] = left.v[ 9] - right.v[ 9];
    result.v[10] = left.v[10] - right.v[10];
    result.v[11] = left.v[11] - right.v[11];
    result.v[12] = left.v[12] - right.v[12];
    result.v[13] = left.v[13] - right.v[13];
    result.v[14] = left.v[14] - right.v[14];
    result.v[15] = left.v[15] - right.v[15];
    return result;
}

Vint16 vshift(const Vint16 &that, const int &val) {
    Vint16 result;
    result.v[15] = that.v[14];
    result.v[14] = that.v[13];
    result.v[13] = that.v[12];
    result.v[12] = that.v[11];
    result.v[11] = that.v[10];
    result.v[10] = that.v[ 9];
    result.v[ 9] = that.v[ 8];
    result.v[ 8] = that.v[ 7];
    result.v[ 7] = that.v[ 6];
    result.v[ 6] = that.v[ 5];
    result.v[ 5] = that.v[ 4];
    result.v[ 4] = that.v[ 3];
    result.v[ 3] = that.v[ 2];
    result.v[ 2] = that.v[ 1];
    result.v[ 1] = that.v[ 0];
    result.v[ 0] = val;
    return result;
}

Vint16 vmax(const Vint16 &left, const Vint16 &right) {
    Vint16 result;
    result.v[ 0] = MAX(left.v[ 0],right.v[ 0]);
    result.v[ 1] = MAX(left.v[ 1],right.v[ 1]);
    result.v[ 2] = MAX(left.v[ 2],right.v[ 2]);
    result.v[ 3] = MAX(left.v[ 3],right.v[ 3]);
    result.v[ 4] = MAX(left.v[ 4],right.v[ 4]);
    result.v[ 5] = MAX(left.v[ 5],right.v[ 5]);
    result.v[ 6] = MAX(left.v[ 6],right.v[ 6]);
    result.v[ 7] = MAX(left.v[ 7],right.v[ 7]);
    result.v[ 8] = MAX(left.v[ 8],right.v[ 8]);
    result.v[ 9] = MAX(left.v[ 9],right.v[ 9]);
    result.v[10] = MAX(left.v[10],right.v[10]);
    result.v[11] = MAX(left.v[11],right.v[11]);
    result.v[12] = MAX(left.v[12],right.v[12]);
    result.v[13] = MAX(left.v[13],right.v[13]);
    result.v[14] = MAX(left.v[14],right.v[14]);
    result.v[15] = MAX(left.v[15],right.v[15]);
    return result;
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
#endif


static int sw(
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
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF;
    printf(" %3d", tbl_pr[0]);
    
    /* first row */
    for (int j=1; j<=s2Len; ++j) {
        tbl_pr[j] = -open -(j-1)*gap;
        del_pr[j] = NEG_INF;
        printf(" %3d", tbl_pr[j]);
    }
    printf("\n");

    /* iter over first sequence */
    for (int i=1; i<=s1Len; ++i) {
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = -open - (i-1)*gap;
        int ins_cr = NEG_INF;
        tbl_pr[0] = Wscore;
        printf(" %3d", Wscore);
        for (int j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + BLOSUM(s1[i-1],s2[j-1]);
            //printf(" NWscore=%d Nscore=%d del=%d ins=%d dig=%d\n",
                    //NWscore, Nscore, del_pr[j], ins_cr, tbl_pr[j]);
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
            printf(" %3d", Wscore);
        }
        printf("\n");
    }

    return tbl_pr[s2Len];
}


/* shift given vector v and insert val */
static inline __m128i vshift16(const __m128i &v, int val)
{
    __m128i ret = _mm_srli_si128(v, 2);
    ret = _mm_insert_epi16(ret, val, 7);
    return ret;
}


static void print_m128i_16(const char *name, const __m128i &m) {
    union {
        __m128i m;
        int16_t v[8];
    } tmp;
    tmp.m = m;
    printf("%s={%3d,%3d,%3d,%3d,%3d,%3d,%3d,%3d}\n", name,
            tmp.v[0], tmp.v[1], tmp.v[2], tmp.v[3],
            tmp.v[4], tmp.v[5], tmp.v[6], tmp.v[7]);
}


static int nw_sse8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr)
{
    int score = 0;
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

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, -open-(i-1)*gap, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(Wscore, 7);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+0)*gap, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(Wscore, 6);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+1)*gap, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(Wscore, 5);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+2)*gap, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(Wscore, 4);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+3)*gap, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(Wscore, 3);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+4)*gap, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(Wscore, 2);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+5)*gap, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(Wscore, 1);
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
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+6)*gap, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = -open-(i+6)*gap;
#if DEBUG
        array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(Wscore, 0);
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
            Wscore = vTbl;
            tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
            del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
            array[(i+0)*(s2Len+14+1) + j + 7] = _mm_extract_epi16(vTbl, 7);
            array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
            array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
            array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
            array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
            array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
            array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
            array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        }
        if (i+0 == s1Len) score = _mm_extract_epi16(vTbl, 7);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+1)*(s2Len+14+1) + j + 6] = _mm_extract_epi16(vTbl, 6);
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+1 == s1Len) score = _mm_extract_epi16(vTbl, 6);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+2)*(s2Len+14+1) + j + 5] = _mm_extract_epi16(vTbl, 5);
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+2 == s1Len) score = _mm_extract_epi16(vTbl, 5);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+3)*(s2Len+14+1) + j + 4] = _mm_extract_epi16(vTbl, 4);
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+3 == s1Len) score = _mm_extract_epi16(vTbl, 4);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+4)*(s2Len+14+1) + j + 3] = _mm_extract_epi16(vTbl, 3);
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+4 == s1Len) score = _mm_extract_epi16(vTbl, 3);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+5)*(s2Len+14+1) + j + 2] = _mm_extract_epi16(vTbl, 2);
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+5 == s1Len) score = _mm_extract_epi16(vTbl, 2);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+6)*(s2Len+14+1) + j + 1] = _mm_extract_epi16(vTbl, 1);
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+6 == s1Len) score = _mm_extract_epi16(vTbl, 1);

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
        Wscore = vTbl;
        tbl_pr[j] = _mm_extract_epi16(vTbl, 0);
        del_pr[j] = _mm_extract_epi16(vDel, 0);
#if DEBUG
        array[(i+7)*(s2Len+14+1) + j + 0] = _mm_extract_epi16(vTbl, 0);
#endif
        if (i+7 == s1Len) score = _mm_extract_epi16(vTbl, 0);
    }

#if DEBUG
    print_array2(array, s1, s1Len, s2, s2Len, 14);
    delete [] array;
#endif
    return score;
}


static int sw_vec2(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    Vint2 Vscore;
    Vint2 Vzero;
    Vint2 Vgap_ext(gap_ext);
    Vint2 Vgap_open(gap_open);
    init_vect(lena+2, nogap, 0);
    init_vect(lena+2, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+2)*(lenb(=%d)+1) = %d\n",
            lena, lenb, (lena+2)*(lenb+1));
    int *array = new int[(lena+2)*(lenb+1)];
    init_vect((lena+2)*(lenb+1), array, -9);
#endif
    for (i=0; i<lenb; i+=2) {
        Vint2 Vmatrix;
        Vint2 Va_gap = -Vgap_open;
        Vint2 Vb_gap(b_gap[0], 0);
        Vint2 Vnogap(nogap[0], 0);
        Vint2 Vlast_nogap = Vzero;
        Vint2 Vprev_nogap = Vzero;
        for (j=0; j<lena+1; ++j) {
            Vb_gap = vshift(Vb_gap, b_gap[j+1]);
            Vnogap = vshift(Vnogap, nogap[j+1]);
            Va_gap = vmax((Vlast_nogap - Vgap_open), (Va_gap - Vgap_ext));
            Vb_gap = vmax((Vnogap - Vgap_open), (Vb_gap - Vgap_ext));
            Vmatrix = Vint2(
                BLOSUM(seqA[j],seqB[i]),
                BLOSUM(seqA[j-1], seqB[i+1])
            );
            Vlast_nogap = vmax((Vprev_nogap + Vmatrix), Vzero);
            Vlast_nogap = vmax(Vlast_nogap, Va_gap);
            Vlast_nogap = vmax(Vlast_nogap, Vb_gap);
            Vprev_nogap = Vnogap;
            Vnogap = Vlast_nogap;
            b_gap[j] = Vb_gap.v[1];
            nogap[j] = Vnogap.v[1];
#if DEBUG
            array[(i+0)*(lena+2) + j+1] = Vlast_nogap.v[0];
            array[(i+1)*(lena+2) + j+0] = Vlast_nogap.v[1];
#endif
            Vscore = vmax(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 2);
    delete [] array;
#endif
    return Vscore.max();
}


static void print_v4(const char *name, const Vint4 &m) {
    printf("%s={%d,%d,%d,%d}\n",
            name, m.v[0], m.v[1], m.v[2], m.v[3]);
}


static int sw_vec4(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    Vint4 Vscore;
    Vint4 Vzero;
    Vint4 Vgap_ext(gap_ext);
    Vint4 Vgap_open(gap_open);
    init_vect(lena+6, nogap, 0);
    init_vect(lena+6, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+6)*(lenb(=%d)+3) = %d\n",
            lena, lenb, (lena+6)*(lenb+3));
    int *array = new int[(lena+6)*(lenb+3)];
    init_vect((lena+6)*(lenb+3), array, -9);
#endif
    for (i=0; i<lenb; i+=4) {
        Vint4 Vmatrix;
        Vint4 Va_gap = -Vgap_open;
        Vint4 Vb_gap(b_gap[2], b_gap[1], b_gap[0], 0);
        Vint4 Vnogap(nogap[2], nogap[1], nogap[0], 0);
        Vint4 Vlast_nogap = Vzero;
        Vint4 Vprev_nogap = Vzero;
        for (j=0; j<lena+3; ++j) {
            Vb_gap = vshift(Vb_gap, b_gap[j+3]);
            Vnogap = vshift(Vnogap, nogap[j+3]);
            Va_gap = vmax((Vlast_nogap - Vgap_open), (Va_gap - Vgap_ext));
            Vb_gap = vmax((Vnogap - Vgap_open), (Vb_gap - Vgap_ext));
            Vmatrix = Vint4(
                BLOSUM(seqA[j-0],seqB[i]),
                BLOSUM(seqA[j-1],seqB[i+1]),
                BLOSUM(seqA[j-2],seqB[i+2]),
                BLOSUM(seqA[j-3],seqB[i+3])
            );
            Vlast_nogap = vmax((Vprev_nogap + Vmatrix), Vzero);
            Vlast_nogap = vmax(Vlast_nogap, Va_gap);
            Vlast_nogap = vmax(Vlast_nogap, Vb_gap);
            Vprev_nogap = Vnogap;
            Vnogap = Vlast_nogap;
            b_gap[j] = Vb_gap.v[3];
            nogap[j] = Vnogap.v[3];
#if DEBUG
            array[(i+0)*(lena+6) + j+3] = Vnogap.v[0];
            array[(i+1)*(lena+6) + j+2] = Vnogap.v[1];
            array[(i+2)*(lena+6) + j+1] = Vnogap.v[2];
            array[(i+3)*(lena+6) + j+0] = Vnogap.v[3];
#endif
            Vscore = vmax(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 6);
    delete [] array;
#endif
    return Vscore.max();
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
    ret.v[0] = ::std::max(vleft.v[0],vright.v[0]);
    ret.v[1] = ::std::max(vleft.v[1],vright.v[1]);
    ret.v[2] = ::std::max(vleft.v[2],vright.v[2]);
    ret.v[3] = ::std::max(vleft.v[3],vright.v[3]);
    return ret.m;
}


/* sse2 doesn't have horizontal max */
static inline int hmax32(const __m128i &m)
{
    int ret;
    v4__m128i tmp;
    tmp.m = m;
    ret = tmp.v[0];
    ret = ::std::max(ret, tmp.v[1]);
    ret = ::std::max(ret, tmp.v[2]);
    ret = ::std::max(ret, tmp.v[3]);
    return ret;
}


static int sw_sse4(
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


static int sw_vec8(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    Vint8 Vscore;
    Vint8 Vzero;
    Vint8 Vgap_ext(gap_ext);
    Vint8 Vgap_open(gap_open);
    init_vect(lena+14, nogap, 0);
    init_vect(lena+14, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+14)*(lenb(=%d)+7) = %d\n",
            lena, lenb, (lena+14)*(lenb+7));
    int *array = new int[(lena+14)*(lenb+7)];
    init_vect((lena+14)*(lenb+7), array, -9);
#endif
    for (i=0; i<lenb; i+=8) {
        Vint8 Vmatrix;
        Vint8 Va_gap = -Vgap_open;
        Vint8 Vb_gap(b_gap[6], b_gap[5], b_gap[4], b_gap[3],
                     b_gap[2], b_gap[1], b_gap[0], 0);
        Vint8 Vnogap(nogap[6], nogap[5], nogap[4], nogap[3],
                     nogap[2], nogap[1], nogap[0], 0);
        Vint8 Vlast_nogap = Vzero;
        Vint8 Vprev_nogap = Vzero;
        for (j=0; j<lena+7; ++j) {
            Vb_gap = vshift(Vb_gap, b_gap[j+7]);
            Vnogap = vshift(Vnogap, nogap[j+7]);
            Va_gap = vmax((Vlast_nogap - Vgap_open), (Va_gap - Vgap_ext));
            Vb_gap = vmax((Vnogap - Vgap_open), (Vb_gap - Vgap_ext));
            Vmatrix = Vint8(
                BLOSUM(seqA[j-0],seqB[i+0]),
                BLOSUM(seqA[j-1],seqB[i+1]),
                BLOSUM(seqA[j-2],seqB[i+2]),
                BLOSUM(seqA[j-3],seqB[i+3]),
                BLOSUM(seqA[j-4],seqB[i+4]),
                BLOSUM(seqA[j-5],seqB[i+5]),
                BLOSUM(seqA[j-6],seqB[i+6]),
                BLOSUM(seqA[j-7],seqB[i+7])
            );
            Vlast_nogap = vmax((Vprev_nogap + Vmatrix), Vzero);
            Vlast_nogap = vmax(Vlast_nogap, Va_gap);
            Vlast_nogap = vmax(Vlast_nogap, Vb_gap);
            Vprev_nogap = Vnogap;
            Vnogap = Vlast_nogap;
            b_gap[j] = Vb_gap.v[7];
            nogap[j] = Vnogap.v[7];
#if DEBUG
            array[(i+0)*(lena+14) + j+7] = Vnogap.v[0];
            array[(i+1)*(lena+14) + j+6] = Vnogap.v[1];
            array[(i+2)*(lena+14) + j+5] = Vnogap.v[2];
            array[(i+3)*(lena+14) + j+4] = Vnogap.v[3];
            array[(i+4)*(lena+14) + j+3] = Vnogap.v[4];
            array[(i+5)*(lena+14) + j+2] = Vnogap.v[5];
            array[(i+6)*(lena+14) + j+1] = Vnogap.v[6];
            array[(i+7)*(lena+14) + j+0] = Vnogap.v[7];
#endif
            Vscore = vmax(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 14);
    delete [] array;
#endif
    return Vscore.max();
}


static int sw_vec16(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    Vint16 Vscore;
    Vint16 Vzero;
    Vint16 Vgap_ext(gap_ext);
    Vint16 Vgap_open(gap_open);
    init_vect(lena+30, nogap, 0);
    init_vect(lena+30, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+30)*(lenb(=%d)+15) = %d\n",
            lena, lenb, (lena+30)*(lenb+15));
    int *array = new int[(lena+30)*(lenb+15)];
    init_vect((lena+30)*(lenb+15), array, -9);
#endif
    for (i=0; i<lenb; i+=16) {
        Vint16 Vmatrix;
        Vint16 Va_gap = -Vgap_open;
        Vint16 Vb_gap(b_gap[14], b_gap[13], b_gap[12], b_gap[11],
                      b_gap[10], b_gap[ 9], b_gap[ 8], b_gap[ 7],
                      b_gap[ 6], b_gap[ 5], b_gap[ 4], b_gap[ 3],
                      b_gap[ 2], b_gap[ 1], b_gap[ 0], 0);
        Vint16 Vnogap(nogap[14], nogap[13], nogap[12], nogap[11],
                      nogap[10], nogap[ 9], nogap[ 8], nogap[ 7],
                      nogap[ 6], nogap[ 5], nogap[ 4], nogap[ 3],
                      nogap[ 2], nogap[ 1], nogap[ 0], 0);
        Vint16 Vlast_nogap = Vzero;
        Vint16 Vprev_nogap = Vzero;
        for (j=0; j<lena+15; ++j) {
            Vb_gap = vshift(Vb_gap, b_gap[j+15]);
            Vnogap = vshift(Vnogap, nogap[j+15]);
            Va_gap = vmax((Vlast_nogap - Vgap_open), (Va_gap - Vgap_ext));
            Vb_gap = vmax((Vnogap - Vgap_open), (Vb_gap - Vgap_ext));
            Vmatrix = Vint16(
                BLOSUM(seqA[j- 0],seqB[i+ 0]),
                BLOSUM(seqA[j- 1],seqB[i+ 1]),
                BLOSUM(seqA[j- 2],seqB[i+ 2]),
                BLOSUM(seqA[j- 3],seqB[i+ 3]),
                BLOSUM(seqA[j- 4],seqB[i+ 4]),
                BLOSUM(seqA[j- 5],seqB[i+ 5]),
                BLOSUM(seqA[j- 6],seqB[i+ 6]),
                BLOSUM(seqA[j- 7],seqB[i+ 7]),
                BLOSUM(seqA[j- 8],seqB[i+ 8]),
                BLOSUM(seqA[j- 9],seqB[i+ 9]),
                BLOSUM(seqA[j-10],seqB[i+10]),
                BLOSUM(seqA[j-11],seqB[i+11]),
                BLOSUM(seqA[j-12],seqB[i+12]),
                BLOSUM(seqA[j-13],seqB[i+13]),
                BLOSUM(seqA[j-14],seqB[i+14]),
                BLOSUM(seqA[j-15],seqB[i+15])
            );
            Vlast_nogap = vmax((Vprev_nogap + Vmatrix), Vzero);
            Vlast_nogap = vmax(Vlast_nogap, Va_gap);
            Vlast_nogap = vmax(Vlast_nogap, Vb_gap);
            Vprev_nogap = Vnogap;
            Vnogap = Vlast_nogap;
            b_gap[j] = Vb_gap.v[15];
            nogap[j] = Vnogap.v[15];
#if DEBUG
            array[(i+ 0)*(lena+30) + j+15] = Vnogap.v[ 0];
            array[(i+ 1)*(lena+30) + j+14] = Vnogap.v[ 1];
            array[(i+ 2)*(lena+30) + j+13] = Vnogap.v[ 2];
            array[(i+ 3)*(lena+30) + j+12] = Vnogap.v[ 3];
            array[(i+ 4)*(lena+30) + j+11] = Vnogap.v[ 4];
            array[(i+ 5)*(lena+30) + j+10] = Vnogap.v[ 5];
            array[(i+ 6)*(lena+30) + j+ 9] = Vnogap.v[ 6];
            array[(i+ 7)*(lena+30) + j+ 8] = Vnogap.v[ 7];
            array[(i+ 8)*(lena+30) + j+ 7] = Vnogap.v[ 8];
            array[(i+ 9)*(lena+30) + j+ 6] = Vnogap.v[ 9];
            array[(i+10)*(lena+30) + j+ 5] = Vnogap.v[10];
            array[(i+11)*(lena+30) + j+ 4] = Vnogap.v[11];
            array[(i+12)*(lena+30) + j+ 3] = Vnogap.v[12];
            array[(i+13)*(lena+30) + j+ 2] = Vnogap.v[13];
            array[(i+14)*(lena+30) + j+ 1] = Vnogap.v[14];
            array[(i+15)*(lena+30) + j+ 0] = Vnogap.v[15];
#endif
            Vscore = vmax(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 30);
    delete [] array;
#endif
    return Vscore.max();
}


typedef union {
    __m128i m;
    int16_t v[8];
} v8__m128i;


static int sw_sse8(
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
            b_gap[j] = _mm_extract_epi16(Vb_gap, 0);
            nogap[j] = _mm_extract_epi16(Vnogap, 0);
#if DEBUG
            array[(i+0)*(lena+14) + j+7] = _mm_extract_epi16(Vnogap, 7);
            array[(i+1)*(lena+14) + j+6] = _mm_extract_epi16(Vnogap, 6);
            array[(i+2)*(lena+14) + j+5] = _mm_extract_epi16(Vnogap, 5);
            array[(i+3)*(lena+14) + j+4] = _mm_extract_epi16(Vnogap, 4);
            array[(i+4)*(lena+14) + j+3] = _mm_extract_epi16(Vnogap, 3);
            array[(i+5)*(lena+14) + j+2] = _mm_extract_epi16(Vnogap, 2);
            array[(i+6)*(lena+14) + j+1] = _mm_extract_epi16(Vnogap, 1);
            array[(i+7)*(lena+14) + j+0] = _mm_extract_epi16(Vnogap, 0);
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


typedef union {
    __m128i m;
    int8_t v[16];
} v16__m128i;


/* shift given vector v and insert val */
static inline __m128i vshift8(const __m128i &v, int val)
{
    v16__m128i tmp;
    tmp.m = _mm_srli_si128(v, 1);
    tmp.v[15] = int8_t(val);
    return tmp.m;
}


static int sw_sse16(
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
    __m128i Vgap_ext = _mm_set1_epi8(gap_ext);
    __m128i Vgap_open = _mm_set1_epi8(gap_open);
    init_vect(lena+30, nogap, 0);
    init_vect(lena+30, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+30)*(lenb(=%d)+15) = %d\n",
            lena, lenb, (lena+30)*(lenb+15));
    int *array = new int[(lena+30)*(lenb+15)];
    init_vect((lena+30)*(lenb+15), array, -9);
#endif
    for (i=0; i<lenb; i+=16) {
        __m128i Vtmp1;
        __m128i Vtmp2;
        __m128i Vmatrix;
        __m128i Va_gap = _mm_sub_epi8(Vzero, Vgap_open);
        __m128i Vb_gap = _mm_set_epi8(b_gap[14], b_gap[13], b_gap[12], b_gap[11],
                                      b_gap[10], b_gap[ 9], b_gap[ 8], b_gap[ 7],
                                      b_gap[ 6], b_gap[ 5], b_gap[ 4], b_gap[ 3],
                                      b_gap[ 2], b_gap[ 1], b_gap[ 0], 0);
        __m128i Vnogap = _mm_set_epi8(nogap[14], nogap[13], nogap[12], nogap[11],
                                      nogap[10], nogap[ 9], nogap[ 8], nogap[ 7],
                                      nogap[ 6], nogap[ 5], nogap[ 4], nogap[ 3],
                                      nogap[ 2], nogap[ 1], nogap[ 0], 0);
        __m128i Vlast_nogap = Vzero;
        __m128i Vprev_nogap = Vzero;
        for (j=0; j<lena+15; ++j) {
            Vb_gap = vshift8(Vb_gap, b_gap[j+15]);
            Vnogap = vshift8(Vnogap, nogap[j+15]);
            Vtmp1 = _mm_sub_epi8(Vlast_nogap, Vgap_open);
            Vtmp2 = _mm_sub_epi8(Va_gap, Vgap_ext);
            Va_gap = _mm_max_epu8(Vtmp1, Vtmp2);
            Vtmp1 = _mm_sub_epi8(Vnogap, Vgap_open);
            Vtmp2 = _mm_sub_epi8(Vb_gap, Vgap_ext);
            Vb_gap = _mm_max_epu8(Vtmp1, Vtmp2);
            Vmatrix = _mm_set_epi8(
                BLOSUM(seqA[j- 0],seqB[i+ 0]),
                BLOSUM(seqA[j- 1],seqB[i+ 1]),
                BLOSUM(seqA[j- 2],seqB[i+ 2]),
                BLOSUM(seqA[j- 3],seqB[i+ 3]),
                BLOSUM(seqA[j- 4],seqB[i+ 4]),
                BLOSUM(seqA[j- 5],seqB[i+ 5]),
                BLOSUM(seqA[j- 6],seqB[i+ 6]),
                BLOSUM(seqA[j- 7],seqB[i+ 7]),
                BLOSUM(seqA[j- 8],seqB[i+ 8]),
                BLOSUM(seqA[j- 9],seqB[i+ 9]),
                BLOSUM(seqA[j-10],seqB[i+10]),
                BLOSUM(seqA[j-11],seqB[i+11]),
                BLOSUM(seqA[j-12],seqB[i+12]),
                BLOSUM(seqA[j-13],seqB[i+13]),
                BLOSUM(seqA[j-14],seqB[i+14]),
                BLOSUM(seqA[j-15],seqB[i+15])
            );
            Vtmp1 = _mm_add_epi8(Vprev_nogap, Vmatrix);
            Vlast_nogap = _mm_max_epu8(Vtmp1, Vzero);
            Vlast_nogap = _mm_max_epu8(Vlast_nogap, Va_gap);
            Vlast_nogap = _mm_max_epu8(Vlast_nogap, Vb_gap);
            Vprev_nogap = Vnogap;
            Vnogap = Vlast_nogap;
            {
                v16__m128i tmp;
                tmp.m = Vb_gap;
                b_gap[j] = tmp.v[15];
                tmp.m = Vnogap;
                nogap[j] = tmp.v[15];
            }
#if DEBUG
            {
                v16__m128i tmp;
                tmp.m = Vnogap;
                array[(i+ 0)*(lena+30) + j+15] = tmp.v[15];
                array[(i+ 1)*(lena+30) + j+14] = tmp.v[14];
                array[(i+ 2)*(lena+30) + j+13] = tmp.v[13];
                array[(i+ 3)*(lena+30) + j+12] = tmp.v[12];
                array[(i+ 4)*(lena+30) + j+11] = tmp.v[11];
                array[(i+ 5)*(lena+30) + j+10] = tmp.v[10];
                array[(i+ 6)*(lena+30) + j+ 9] = tmp.v[ 9];
                array[(i+ 7)*(lena+30) + j+ 8] = tmp.v[ 8];
                array[(i+ 8)*(lena+30) + j+ 7] = tmp.v[ 7];
                array[(i+ 9)*(lena+30) + j+ 6] = tmp.v[ 6];
                array[(i+10)*(lena+30) + j+ 5] = tmp.v[ 5];
                array[(i+11)*(lena+30) + j+ 4] = tmp.v[ 4];
                array[(i+12)*(lena+30) + j+ 3] = tmp.v[ 3];
                array[(i+13)*(lena+30) + j+ 2] = tmp.v[ 2];
                array[(i+14)*(lena+30) + j+ 1] = tmp.v[ 1];
                array[(i+15)*(lena+30) + j+ 0] = tmp.v[ 0];
            }
#endif
            Vscore = _mm_max_epu8(Vscore, Vlast_nogap);
        }
    }
#if DEBUG
    print_array(array, seqA, lena, seqB, lenb, 30);
    delete [] array;
#endif
    max16(score, Vscore);
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
    int score;
    unsigned long long timer;
#if DEBUG
    size_t limit = 1;
#else
    size_t limit = 1000;
#endif
    size_t i;
    
    timer_init();
    ::std::cout << timer_name() << " timer" << ::std::endl;

    ::std::cout << "alg\t\tscore\ttime" << ::std::endl;

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw orig\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "nw orig\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "nw_sse8\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_vec2(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_vec2\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_vec4(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_vec4\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_vec8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_vec8\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_vec16(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_vec16\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_sse4(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_sse4\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_sse8(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_sse8\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

    /* doesn't work yet, needs to be all unsigned */
#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_sse16(&seqA[16], lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    }
    timer = timer_end(timer);
    ::std::cout << "sw_sse16\t\t" << score << "\t" << timer/limit << ::std::endl;
#endif

    delete [] nogap;
    delete [] b_gap;

    return 0;
}
