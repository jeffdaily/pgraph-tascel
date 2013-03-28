#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "bucket.h"
#include "cfg.h"
#include "elib.h"
#include "lib.h"
#include "loadseq.h"
#include "pairs.h"
#include "search.h"
#include "stree.h"
#include "uFind.h"

#define MAX_FILENAME_LEN 200

static char gSeqFile[MAX_FILENAME_LEN];
static char gCfgFile[MAX_FILENAME_LEN];
static size_t gN;

static int optCK(int argc, char **argv);

int main(int argc, char *argv[]){
    sequence_t *seqs = NULL;
    size_t maxSeqLen;
    size_t NN;
    bucket_t *bkt = NULL;     /* counters for *sf */
    int bktSize;
    suffix_t *sf = NULL;
    int sfSize;
    ufind_t *uSet = NULL; 
    int singletons = 2;
    param_t param;

    setbuf(stdout, (char *)0);

    /* reading in pamameters */
    if(optCK(argc, argv) != 3){
        printf("Usage : %s -f {fasta} -n {#seqs} -c {cfg file}\n", argv[0]);
        exit(0);
    }

    /* read in configurations */
    param.window_size = get_config_val(gCfgFile, "SlideWindowSize");
    param.exact_match_len = get_config_val(gCfgFile, "ExactMatchLen");
    param.AOL = get_config_val(gCfgFile, "AlignOverLongerSeq");
    param.SIM = get_config_val(gCfgFile, "MatchSimilarity");
    param.OS = get_config_val(gCfgFile, "OptimalScoreOverSelfScore");
    
    init_map(SIGMA);

    printf("----------------------------------------------\n");
    printf("%-15s: %.32s\n", "fasta seq", gSeqFile);
    printf("%-15s: %.32s\n", "config file", gCfgFile);
    printf("%-15s: %d\n", "slide size", param.window_size);
    printf("%-15s: %d\n", "exactMatch len", param.exact_match_len);
    printf("----------------------------------------------\n");

    seqs = emalloc(gN*(sizeof *seqs));

    time_t t1, t2;
    /*------------*
     * load seqs 
     *------------*/
    (void) time(&t1);
    load_all_sequences(gSeqFile, gN, seqs, &NN, &maxSeqLen);
    (void) time(&t2);
    printf("<%zu> amino acids are loaded in <%lld> secs\n",
            NN, (long long)(t2-t1));



    /* allocate mem for *bkt, *sf */
    bktSize = power(SIGMA, param.window_size); /* bucket size */
    bkt = emalloc(bktSize*(sizeof *bkt));

    sfSize = NN - gN*param.window_size; /* NN includes '$' */
    sf = emalloc(sfSize*(sizeof *sf));


    /*-------------------*
     * build the bucket 
     *-------------------*/ 
    (void) time(&t1);
    build_buckets(seqs, gN, bkt, bktSize, sf, sfSize, param.window_size); 
    (void) time(&t2);
    printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));

    
    (void) time(&t1);
    cntSort4Bkt(bkt, bktSize);
    //qsort(bkt, bktSize, sizeof *bkt, bktCmp);
    (void) time(&t2);
    printf("Bucketing sorted in <%lld> secs\n", (long long)(t2-t1));


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
    printf("Tree constructed and processed in <%lld> secs\n",
            (long long)(t2-t1));
    print_pairs();

    
    disp_all_clusters(uSet, gN, &singletons, ".", seqs); 

    /*-----------*
     * free mem 
     *-----------*/
    free(sf);
    free(bkt);

    free_sequences(seqs, gN);

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
                gN = atol(optarg);
			    cnt++;
			    break;
		}
	}
	return cnt;
}
