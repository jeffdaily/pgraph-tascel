#include <stdlib.h>
#include <stdio.h>

#if HAVE_GMP
#include <gmp.h>
#endif

#include "combinations.h"


unsigned long binomial_coefficient(unsigned long n, unsigned long k)
{
#if HAVE_GMP
    mpz_t mpz_answer;
    unsigned long ui_answer;

    if (n<k) {
        return 0;
    }
    if (n == k) {
        return 1;
    }

    mpz_init(mpz_answer);
    mpz_bin_uiui(mpz_answer, n, k);
    ui_answer = mpz_get_ui(mpz_answer);
    mpz_clear(mpz_answer);

    return ui_answer;
#else
    /* from http://blog.plover.com/math/choose.html */
    unsigned long r = 1;
    unsigned long d;
    if (k > n) return 0;
    for (d=1; d <= k; d++) {
        r *= n--;
        r /= d;
    }
    return r;
#endif
}


void k_combination(unsigned long pos, unsigned long k, unsigned long *result)
{
    unsigned long n;
    unsigned long i;
    unsigned long bc;
    unsigned long bc_previous;

    for (i=k; i>0; --i) {
        if (0 == pos) {
            result[i-1] = i-1;
            continue;
        }
        n=i-1;
        bc_previous = 0;
        bc = binomial_coefficient(n,i);
        while (bc <= pos) {
            bc_previous = bc;
            ++n;
            bc = binomial_coefficient(n,i);
        }
        pos -= bc_previous;
        result[i-1] = n-1;
    }
}


void init_combination(unsigned long k, unsigned long *combination)
{
    size_t i;

    for (i=0; i<k; ++i) {
        combination[i] = i;
    }
}


void next_combination(unsigned long k, unsigned long *combination)
{
    size_t i;

    /* as we find position to increment, reset values */
    for (i=0; i<k-1; ++i) {
        if (combination[i]+1 < combination[i+1]) {
            break;
        }
        else {
            combination[i] = i;
        }
    }
    /* increment it */
    ++combination[i];
}
