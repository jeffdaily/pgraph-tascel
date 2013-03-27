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


#ifdef __cplusplus
extern "C" {
#endif

void loadAllSeqs(char *fileName, int nSeqs, SEQ *seqs, long *nChars);
void printSeq(SEQ *seqs, int index);
int freeSeqs(SEQ *seqs, int nSeqs);

#ifdef __cplusplus
}
#endif

#endif /* end of loadseq.h */
