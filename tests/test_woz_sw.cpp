#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blosum/blosum62.h"

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

#define BLOSUM(ch1, ch2) (matrix[MAP_BLOSUM_[(ch1)]][MAP_BLOSUM_[(ch2)]])
#define MAX(a,b) ((a)>(b)?(a):(b))
#define DEBUG 1

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


static int sw(
        const char *seqA, const int lena,
        const char *seqB, const int lenb,
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
    printf(" ");
    for (j=0; j<lena; ++j) {
        printf("  %c", seqA[j]);
    }
    printf("\n");
#endif
    for (i=0; i<lenb; ++i) {
        int a_gap;
        last_nogap = prev_nogap = 0;
        a_gap = -gap_open;
#if DEBUG
        printf("%c", seqB[i]);
#endif
        for (j=0; j<lena; ++j) {
            a_gap = MAX((last_nogap - gap_open), (a_gap - gap_ext));
            b_gap[j] = MAX((nogap[j] - gap_open), (b_gap[j] - gap_ext));
            last_nogap = MAX((prev_nogap + BLOSUM(seqA[j],seqB[i])), 0);
            last_nogap = MAX(last_nogap, a_gap);
            last_nogap = MAX(last_nogap, b_gap[j]);
            prev_nogap = nogap[j];
            nogap[j] = last_nogap;
            printf(" %2d", last_nogap);
            score = MAX(score, last_nogap);
        }
        printf("\n");
    }
    return score;
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
    int padb = pad/2 + 1;

    printf("   $");
    for (j=1; j<pad/2; ++j) {
        printf("  $");
    }
    for (j=0; j<lena; ++j) {
        printf("  %c", seqA[j]);
    }
    for (j=0; j<pad/2; ++j) {
        printf("  $");
    }
    printf("\n");
    for (i=0; i<lenb; ++i) {
        printf("%c", seqB[i]);
        for (j=0; j<lena+pad; ++j) {
            printf(" %2d", array[i*(lena+pad) + j]);
        }
        printf("\n");
    }
    for (i=lenb; i<(lenb+lenb%padb); ++i) {
        printf("$");
        for (j=0; j<lena+pad; ++j) {
            printf(" %2d", array[i*(lena+pad) + j]);
        }
        printf("\n");
    }
}
#endif


static int sw_vec2(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    int score;
    Vint2 Vscore;
    Vint2 Vzero;
    Vint2 Vgap_ext(gap_ext);
    Vint2 Vgap_open(gap_open);
    init_vect(lena+2, nogap, 0);
    init_vect(lena+2, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+2)*lenb(=%d) = %d\n",
            lena, lenb, (lena+2)*lenb);
    int *array = new int[(lena+2)*lenb];
    init_vect((lena+2)*lenb, array, -9);
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
            if (i==0 && j==0) {
                printf("Vb_gap={%d,%d}\n", Vb_gap.v[0], Vb_gap.v[1]);
            }
            Vnogap = vshift(Vnogap, nogap[j+1]);
            if (i==0 && j==0) {
                printf("Vnogap={%d,%d}\n", Vnogap.v[0], Vnogap.v[1]);
            }
            Va_gap = vmax((Vlast_nogap - Vgap_open), (Va_gap - Vgap_ext));
            if (i==0 && j==0) {
                printf("Va_gap={%d,%d}\n", Va_gap.v[0], Va_gap.v[1]);
            }
            Vb_gap = vmax((Vnogap - Vgap_open), (Vb_gap - Vgap_ext));
            if (i==0 && j==0) {
                printf("Vb_gap={%d,%d}\n", Vb_gap.v[0], Vb_gap.v[1]);
            }
            Vmatrix = Vint2(
                BLOSUM(seqA[j],seqB[i]),
                BLOSUM(seqA[j-1], seqB[i+1])
            );
            if (i==0 && j==0) {
                printf("Vmatrix={%d,%d}\n", Vmatrix.v[0], Vmatrix.v[1]);
            }
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
            if (i==0) {
                printf("Vlast_nogap={%d,%d}\n", Vlast_nogap.v[0], Vlast_nogap.v[1]);
            }
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


static int sw_vec4(
        const char * const restrict seqA, const int lena,
        const char * const restrict seqB, const int lenb,
        const int gap_open, const int gap_ext,
        const int matrix[24][24],
        int * const restrict nogap, int * const restrict b_gap)
{
    int i;
    int j;
    int score;
    Vint4 Vscore;
    Vint4 Vzero;
    Vint4 Vgap_ext(gap_ext);
    Vint4 Vgap_open(gap_open);
    init_vect(lena+6, nogap, 0);
    init_vect(lena+6, b_gap, -gap_open);
#if DEBUG
    printf("array length (lena(=%d)+6)*(lenb(=%d)+lenb%4) = %d\n",
            lena, lenb, (lena+6)*(lenb+lenb%4));
    int *array = new int[(lena+6)*(lenb+lenb%4)];
    init_vect((lena+6)*(lenb+lenb%4), array, -9);
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


int main(int argc, char **argv)
{
    const char *seqA = "SLPSMRADSFTKELMEKISS";
    const char *seqB = "MTNKICIYAISKNEEKFV";
    int lena = strlen(seqA);
    int lenb = strlen(seqB);
    int lenmax = MAX(lena,lenb)+16;
    int *nogap = new int[lenmax];
    int *b_gap = new int[lenmax];
    int score;
    
    score = sw(seqA, lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    printf("%s\n%s\n%d\n", seqA, seqB, score);
    score = sw_vec2(seqA, lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    printf("%s\n%s\n%d\n", seqA, seqB, score);
    score = sw_vec4(seqA, lena, seqB, lenb, 10, 1, blosum62, nogap, b_gap);
    printf("%s\n%s\n%d\n", seqA, seqB, score);
    
    delete [] nogap;
    delete [] b_gap;

    return 0;
}
