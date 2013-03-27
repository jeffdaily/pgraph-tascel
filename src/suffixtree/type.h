#ifndef TYPE_H_
#define TYPE_H_

#define DOLLAR 'U' /* dollar delimiter */
#define BEGIN 'O'  /* it is oh, not zero */
#define SIGMA 26   /* size of alphabet, '*'->'J' */

typedef unsigned int uint;
typedef unsigned long ulong;
enum {NO = 0, YES = 1};
enum {FALSE = 0, TRUE = 1};


/**
 * seq info for fasta file
 */
typedef struct seq{
    char *gid;  /* gid of fasta file */
    char *str;  /* actual string of fasta file */
    int strLen; /* string length */
}SEQ;

/**
 * lookup table, used to represent suffix of 
 * every sequences
 */
typedef struct suff{
    int sid;           /* string id */
    int pid;           /* position id */
    struct suff *next; /* ptr to next suff */
}SUFFIX;

/**
 * bucket structure
 */
typedef struct bkt{
    SUFFIX *bktList;
    int bktCnt;
}BKT;

/**
 * suffix tree node
 */
typedef struct stnode{
    int depth;                   /* depth since the root, not including the initial size k */
    int rLeaf;                  /* right most leaf index */
    struct suff *lset[SIGMA];    /* subtree's nodes branched according to left characters */
}STNODE;

/**
 * CELL for dynamic alignment
 */
typedef struct cell{
    int score;           /* alignment score */
    int ndig;            /* #(matches) */
    int alen;            /* alignment length */
}CELL;

/**
 * parameters packed in struct
 */
typedef struct param{
    int k;               /* slide window size */
    int AOL;             /* AlignOverLongerSeq */
    int SIM;             /* MatchSimilarity */
    int OS;              /* OptimalScoreOverSelfScore */
    int exactMatchLen;   /* exact Match Length cutoff */
}PARAM;

/**
 * statistics packed in struct 
 */
typedef struct work{
    int nAlign;
    int nPairs;
}WORK;

#endif /* end of type.h */
