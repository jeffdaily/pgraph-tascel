/**
 * @file combinations.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _COMBINATIONS_H_
#define _COMBINATIONS_H_

#ifdef __cplusplus
extern "C" {
#endif

unsigned long binomial_coefficient(unsigned long n, unsigned long k);
void k_combination(unsigned long pos, unsigned long k, unsigned long *result);
unsigned long k_combination_inv(unsigned long k, unsigned long *result);
void k_combination2(unsigned long pos, unsigned long *result);
unsigned long k_combination2_inv(unsigned long *result);
void init_combination(unsigned long k, unsigned long *combination);
void next_combination(unsigned long k, unsigned long *combination);
void inc_combination(unsigned long inc, unsigned long k, unsigned long *combination);
void inc_combination2(unsigned long inc, unsigned long *combination);

#ifdef __cplusplus
}
#endif

#endif /* _COMBINATIONS_H_ */
