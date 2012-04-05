#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

#include "combinations.h"

int main(int argc, char **argv)
{
    unsigned long i;
    unsigned long pos = 72;
    unsigned long k = 5;
    unsigned long *answer = NULL;
    unsigned long *series = NULL;

    answer = (unsigned long*)malloc(sizeof(unsigned long) * k);
    k_combination(pos,k,answer);
    printf("%lu={%lu", pos, answer[0]);
    for (i=1; i<k; ++i) {
        printf(",%lu", answer[i]);
    }
    printf("}\n");
    free(answer);

    answer = (unsigned long*)malloc(sizeof(unsigned long) * k);
    series = (unsigned long*)malloc(sizeof(unsigned long) * k);
    for (i=0; i<k; ++i) {
        series[i] = i;
    }
    int counter = 0;
    for (i=0; i<126; ++i) {
        unsigned int j;
        k_combination(counter, k, answer);
        cout << setw(4) << counter++ << ": " << series[0];
        for (j=1; j<k; ++j) {
            cout << "," << series[j];
        }
        cout << " " << answer[0];
        for (j=1; j<k; ++j) {
            cout << "," << answer[j];
        }
        cout << endl;
        /* verify each element is the same */
        for (j=0; j<k; ++j) {
            if (answer[j] != series[j]) {
                cout << "FAILURE" << endl;
                return EXIT_FAILURE;
            }
        }
        next_combination(k, series);
    }

    pos = 204799679999;
    k = 2;
    answer = (unsigned long*)malloc(sizeof(unsigned long) * k);
    k_combination(pos,k,answer);
    printf("%lu={%lu", pos, answer[0]);
    for (i=1; i<k; ++i) {
        printf(",%lu", answer[i]);
    }
    printf("}\n");
    free(answer);

    return 0;
}
