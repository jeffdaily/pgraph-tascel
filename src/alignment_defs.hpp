/**
 * @file alignment.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Routines for dynamic programming of sequence alignment. These routines
 * strive to be as performance-optimal as possible as well as only using
 * standard C data types rather than higher-level C structs. Higher-level C
 * structs would likely simplify the otherwise long argument lists to many of
 * these routines.
 */
#ifndef _PGRAPH_ALIGNMENT_DEFS_H_
#define _PGRAPH_ALIGNMENT_DEFS_H_

/** alphabet for BLOSUM in the order in which the matrix is defined */
static char ALPHABET_BLOSUM[] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        'B', 'Z', 'X', 'J', 'O', 'U' };


/** mapping from BLOSUM alphabet to BLOSUM index; use as BLOSUM[map[ch-'A']] */
static int MAP_BLOSUM[] = {
          0,  20,   4,   3,   6,  13,   7,   8,   9,  23,
         11,  10,  12,   2,  24,  14,   5,   1,  15,  16,
         25,  19,  17,  22,  18,  21 };


static int ALPHABET_DNA[] = { 'A', 'C', 'G', 'T' };


static int MAP_DNA[] = {
          0,  -1,   1,  -1,  -1,  -1,   2,  -1,  -1,  -1,
         -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   3,
         -1,  -1,  -1,  -1,  -1,  -1 };


static int SUB_DNA[4][4] = {
        /*       A   C   G   T */
        /* A */{ 1,  0,  0,  0},
        /* C */{ 0,  1,  0,  0},
        /* G */{ 0,  0,  1,  0},
        /* T */{ 0,  0,  0,  1}
};

#endif /* _PGRAPH_ALIGNMENT_DEFS_H_ */
