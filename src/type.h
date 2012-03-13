/*
 * $Rev: 77 $ 
 * $Date: 2011-05-11 16:03:49 -0700 (Wed, 11 May 2011) $ 
 * $Author: andy.cj.wu@gmail.com $
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * ----------------------------------------------------------------
 *
 */

#ifndef TYPE_H_
#define TYPE_H_

#include <mpi.h>

#define COLOR_PD 0
#define COLOR_SM 1
#define COLOR_MS 2
#define COLOR_CS 3

#define MSG_PM_TAG 1
#define MSG_MC_TAG 2
#define MSG_CM_TAG 3
#define MSG_SM_TAG 4
#define MSG_PS_TAG 5
#define MSG_SP_TAG 6
#define MSG_CR_TAG 7   /* seq request tag */
#define MSG_CS_TAG 8   /* seq string tag */


/* tag to differentiate msg flowing in the system */
#define TAG_P 'P'  /* pair workload from producer */
#define TAG_C 'C'  /* Job request from consumer */
#define TAG_E 'E'  /* ending of pairs */
#define TAG_S 'S'  /* stop signal for program */
#define TAG_T 'T'  /* ending signal of subgroup */
#define TAG_F 'F'  /* end signal for tree file */
#define TAG_K 'K'  /* final end singla */


/* for consumer request */
#define R_SIZE 3
#define R_HALF 2
#define R_QUAT 1
#define R_NONE 0


#define DOLLAR 'U' /* dollar delimiter */
#define BEGIN 'O'  /* it is oh, not zero */
#define SIGMA 26   /* size of alphabet, '*'->'J' */


#define MAX_FILENAME_LEN 200
#define FILE_STOP -1


/* seq. str. stat */
#define SEQ_N 0  /* seq. str is empty */
#define SEQ_S 1  /* statically cached locally */
#define SEQ_R 2  /* req. sent out already */
#define SEQ_F 3  /* seq. str fetched locally */

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long u64;
enum {NO = 0, YES = 1};
enum {FALSE = 0, TRUE = 1};


/* NOTE: if you are changing this struct,
   make sure to change the dtype.h accordingly
 */
typedef struct msg {
    char tag;    /* 'P', 'C', 'E' or 'A' */
    int id1;     /* s1 & rank */
    int id2;     /* s2 & status */
}MSG;


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
 * suffix tree node
 */
typedef struct stnode{
    int depth;                   /* depth since the root, not including the initial size k */
    int rLeaf;                  /* right most leaf index */
    struct suff *lset[SIGMA];    /* subtree's nodes branched according to left characters */
}STNODE;


/**
 * seq info for fasta file
 */
typedef struct seq{
    
    /* isReady TAG
     * ----------------------------------------------------------------
       0: not ready       #define SEQ_N 0  - seq. str is empty
       1: precached       #define SEQ_S 1  - statically cached locally
       2: req. sent out   #define SEQ_R 2  - req. sent out already
       3: seq. str ready  #define SEQ_F 3  - seq. str fetched locally
     */
    char stat;      /* seq string is ready for not? */
    char *str;      /* actual string of fasta file */
    int strLen;     /* string length */
    int cnt;        /* counter for seq usage */
}SEQ;

#endif /* end of type.h */
