#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

#include "combinations.h"

static void print_diff(
        unsigned long k, unsigned long *left, unsigned long *right)
{
    unsigned i;
    cout << (long long)left[0] - (long long)right[0];
    for (i=1; i<k; ++i) {
        cout << "," << (long long)left[i] - (long long)right[i];
    }
}

static void print_combo(unsigned long k, unsigned long *combo)
{
    unsigned long i;
    cout << combo[0];
    for (i=1; i<k; ++i) {
        cout << "," << combo[i];
    }
}

int test(unsigned long k, unsigned long step, unsigned long limit)
{
    unsigned long i = 0;
    unsigned long pos = 0;
    unsigned long *comb1 = NULL;
    unsigned long *comb2 = NULL;
    unsigned long *comb2_prev = NULL;
    unsigned long *comb3 = NULL;

    /* test k_combination, next_combination, inc_combination functions */
    comb1 = (unsigned long*)malloc(sizeof(unsigned long) * k);
    comb2 = (unsigned long*)malloc(sizeof(unsigned long) * k);
    comb2_prev = (unsigned long*)malloc(sizeof(unsigned long) * k);
    comb3 = (unsigned long*)malloc(sizeof(unsigned long) * k);
    init_combination(k, comb2);
    init_combination(k, comb2_prev);
    init_combination(k, comb3);
    for (i=0; i<limit; ++i) {
        unsigned int j;
        if (i%step == 0 && i!=0) {
            inc_combination(step, k, comb3);
        }
        k_combination(i, k, comb1);
        cout << setw(4) << i << ": ";
        print_combo(k, comb2);
        cout << " ";
        print_combo(k, comb1);
        cout << " ";
        print_combo(k, comb3);
        cout << " ";
        print_diff(k, comb2, comb2_prev);
        cout << endl;
        /* verify each element is the same */
        for (j=0; j<k; ++j) {
            if (comb1[j] != comb2[j]) {
                cout << "FAILURE" << endl;
                return EXIT_FAILURE;
            }
            if (i%step == 0 && i!=0 && comb2[j] != comb3[j]) {
                cout << "FAILURE" << endl;
                return EXIT_FAILURE;
            }
        }
        for (j=0; j<k; ++j) {
            comb2_prev[j] = comb2[j];
        }
        next_combination(k, comb2);
    }
    free(comb1);
    free(comb2);
    free(comb2_prev);
    free(comb3);

    /* test a really huge k_combination */
    pos = 204799679999;
    comb1 = (unsigned long*)malloc(sizeof(unsigned long) * k);
    k_combination(pos,k,comb1);
    printf("%lu={%lu", pos, comb1[0]);
    for (i=1; i<k; ++i) {
        printf(",%lu", comb1[i]);
    }
    printf("}\n");
    free(comb1);

    return 0;
}

int main(int argc, char **argv)
{
    int retcode;

    retcode = test(5, 3, 127);

    if (retcode != EXIT_FAILURE) {
        retcode = test(5, 11, 127);
    }

    if (retcode != EXIT_FAILURE) {
        retcode = test(5, 111, 127);
    }

    if (retcode != EXIT_FAILURE) {
        retcode = test(2, 3, 127);
    }

    if (retcode != EXIT_FAILURE) {
        retcode = test(2, 11, 127);
    }

    if (retcode != EXIT_FAILURE) {
        retcode = test(2, 111, 12700);
    }

    return retcode;
}
