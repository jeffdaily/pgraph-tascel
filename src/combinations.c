#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#if HAVE_GMP
#include <gmp.h>
#endif

#include <math.h>

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

    if (2 == k) {
        k_combination2(pos, result);
    }
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


void k_combination2(unsigned long pos, unsigned long *result)
{
    double s;
    double i = floor(sqrt(2.0*pos)) - 1.0;
    if (i <= 1) {
        i = 1.0;
    }
    s = i*(i-1.0)/2.0;
    while (pos-s >= i) {
        s += i;
        i += 1;
    }
    result[0] = pos-s;
    result[1] = i;
}


void init_combination(unsigned long k, unsigned long *combination)
{
    unsigned long i;

    for (i=0; i<k; ++i) {
        combination[i] = i;
    }
}


void next_combination(unsigned long k, unsigned long *combination)
{
    unsigned long i;

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


void inc_combination(
        unsigned long inc, unsigned long k, unsigned long *combination)
{
    unsigned long pos = 0;
    unsigned long i = 0;

    /* special case 2 combination for faster compute */
    if (2 == k) {
        inc_combination2(inc, combination);
    }
    else {
        /* find starting position, reset values */
        for (pos=0; pos<k-1; ++pos) {
            if (combination[pos]+1 < combination[pos+1]) {
                break;
            }
            else {
                combination[pos] = pos;
            }
        }

        while (i < inc) {
            if (pos == k-1 || combination[pos]+1 < combination[pos+1]) {
                /* move all left */
                while (pos>0 && i<inc) {
                    ++i;
                    ++combination[pos];
                    --pos;
                }
                if (i<inc) {
                    ++i;
                    ++combination[pos];
                }
            }
            else {
                while (pos<k-1 && combination[pos]+1 == combination[pos+1]) {
                    combination[pos] = pos;
                    ++pos;
                }
            }
        }
    }
}


#if 1
void inc_combination2(unsigned long inc, unsigned long *combination)
{
    unsigned long i = 0;

    //assert(combination[0] < combination[1]);
    while (i < inc) {
        unsigned long diff = combination[1] - combination[0];
        if ((inc-i) < diff) { 
            combination[0] += inc-i;
            break;
        }
        else {
            i += diff;
            combination[0] = 0;
            ++combination[1];
        }
    }
}
#else
/* this implementation was easier to understand but turned out to be slower */
void inc_combination2(unsigned long inc, unsigned long *combination)
{
    unsigned long diff = combination[1] - combination[0];
    if (inc >= diff) {
        combination[0] = 0;
        ++combination[1];
        inc -= diff;
    }
    while (inc >= combination[1]) {
        inc -= combination[1];
        ++combination[1];
    }
    combination[0] += inc;
}
#endif
