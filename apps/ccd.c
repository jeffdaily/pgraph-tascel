#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "type.h"
#include "loadseq.h"
#include "type.h"
#include "lib.h"
#include "elib.h"
#include "bucket.h"
#include "stree.h"
#include "cfg.h"
#include "search.h"
#include "uFind.h"

#define MAX_FILENAME_LEN 200

static char gSeqFile[MAX_FILENAME_LEN];
static char gCfgFile[MAX_FILENAME_LEN];
static int gN;

#pragma mta parallel off
int optCK(int argc, char **argv);

int main(int argc, char *argv[]){
    SEQ *seqs = NULL;
    STNODE *stNodes = NULL;
    int nStNodes;
    int maxSeqLen;
    int NN;
    BKT *bkt = NULL;     /* counters for *sf */
    int bktSize;
    SUFFIX *sf = NULL;
    int sfSize;
    int exactMatchLen;   /* exact match length cutoff */
    UF *uSet = NULL; 
    int singletons = 2;
    PARAM param;

    setbuf(stdout, (char *)0);

    /* reading in pamameters */
    if(optCK(argc, argv) != 3){
        printf("Usage : %s -f {fasta} -n {#seqs} -c {cfg file}\n", argv[0]);
        exit(0);
    }

    /* read in configurations */
    param.k = getCfgVal(gCfgFile, "SlideWindowSize");   /* window size */
    param.exactMatchLen = getCfgVal(gCfgFile, "ExactMatchLen");      /* exactMatchLen */ 
    param.AOL = getCfgVal(gCfgFile, "AlignOverLongerSeq");
    param.SIM = getCfgVal(gCfgFile, "MatchSimilarity");
    param.OS = getCfgVal(gCfgFile, "OptimalScoreOverSelfScore");
    
    initMap(SIGMA);

    printf("----------------------------------------------\n");
    printf("%-15s: %.32s\n", "fasta seq", gSeqFile);
    printf("%-15s: %.32s\n", "config file", gCfgFile);
    printf("%-15s: %d\n", "slide size", param.k);
    printf("%-15s: %d\n", "exactMatch len", param.exactMatchLen);
    printf("----------------------------------------------\n");

    seqs = emalloc(gN*(sizeof *seqs));

    time_t t1, t2;
    /*------------*
     * load seqs 
     *------------*/
    (void) time(&t1);
    loadAllSeqs(gSeqFile, gN, seqs, &NN, &maxSeqLen);
    (void) time(&t2);
    printf("<%ld> amino acids are loaded in <%d> secs\n", NN, (int)t2-t1);



    /* allocate mem for *bkt, *sf */
    bktSize = power(SIGMA, param.k); /* bucket size */
    bkt = emalloc(bktSize*(sizeof *bkt));

    sfSize = NN - gN*param.k; /* NN includes '$' */
    sf = emalloc(sfSize*(sizeof *sf));


    /*-------------------*
     * build the bucket 
     *-------------------*/ 
    (void) time(&t1);
    buildBkt(seqs, gN, bkt, bktSize, sf, sfSize, param.k); 
    (void) time(&t2);
    printf("Bucketing finished in <%d> secs\n", (int)t2-t1);

    
    (void) time(&t1);
    cntSort4Bkt(bkt, bktSize);
    //qsort(bkt, bktSize, sizeof *bkt, bktCmp);
    (void) time(&t2);
    printf("Bucketing sorted in <%d> secs\n", (int)t2-t1);


    #ifdef DEBUG
    for(int i = 0; i < 100; i++){
        printf("bktSize[%d] - %d\n", i, bkt[i].bktCnt);
    }
    exit(0);
    #endif

    uSet = init_union(gN);

    /*---------------------------------------*
     * suffix tree construction & processing  
     *---------------------------------------*/
    (void) time(&t1);
    buildForest(bkt, bktSize, seqs, gN, maxSeqLen, uSet, &param);
    (void) time(&t2);
    printf("Tree constructed and processed in <%d> secs\n", (int)t2-t1);
    printPairs();

    
    disp_all_clusters(uSet, gN, &singletons, ".", seqs); 

    /*-----------*
     * free mem 
     *-----------*/
    free(sf);
    free(bkt);

    #ifndef CRAY_XMT
    freeSeqs(seqs, gN);
    #endif 

    free(seqs);
    free_union(uSet);
    return EXIT_SUCCESS; 
}

int optCK(int argc, char **argv){
	int option;
	int cnt = 0;
	
	while(-1 != (option = getopt(argc, argv, "f:n:c:"))){
		switch (option){
		    case '?':
                printf("Usage : %s -f {fasta} -n {#seqs} -c {cfg file}\n", argv[0]);
			    exit(-1);
		    case 'f':
			    strcpy(gSeqFile, optarg);
			    cnt++;
			    break;
            case 'c':
                strcpy(gCfgFile, optarg);
                cnt++;
                break;
		    case 'n':
                gN = atoi(optarg);
			    cnt++;
			    break;
		}
	}
	return cnt;
}
