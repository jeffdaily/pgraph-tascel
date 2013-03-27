#ifndef LOAD_SEQ_H_
#define LOAD_SEQ_H_

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "type.h"
#include "elib.h"
#include "wlib.h"

#define MAX_LINE_LEN 100000
#define FASTA_FLAG '>'

#ifdef CRAY_XMT

#ifdef __cplusplus
extern "C" {
#endif

void loadAllSeqs(char *fileName, int nSeqs, SEQ *seqs, int *nChars, int *maxSeqLen);
void printSeq(SEQ *seqs, int index);
int freeSeqs(SEQ *seqs, int nSeqs);

#ifdef __cplusplus
}

#endif  /* end of __c++ */
#endif /* end of CRAY_XMT */

#ifndef CRAY_XMT
void loadAllSeqs(char *fileName, int nSeqs, SEQ *seqs, int *nChars, int *maxSeqLen);
void printSeq(SEQ *seqs, int index);
int freeSeqs(SEQ *seqs, int nSeqs);
#endif 

#endif /* end of loadseq.h */
