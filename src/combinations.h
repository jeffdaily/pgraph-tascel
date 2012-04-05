#ifndef _COMBINATIONS_H_
#define _COMBINATIONS_H_

#ifdef __cplusplus
extern "C" {
#endif 

unsigned long binary_coefficient(unsigned long n, unsigned long k);
void k_combination(unsigned long pos, unsigned long k, unsigned long *result);
void next_combination(unsigned long k, unsigned long *combination);

#ifdef __cplusplus
}
#endif 

#endif /* _COMBINATIONS_H_ */
