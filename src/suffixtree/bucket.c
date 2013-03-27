#include "bucket.h"

static int sfIndex = 0;

/* -----------------------------------------* 
 * Bkt here is essentially counters. Init 
 * them into zero.
 *
 * @param bkt - headers for *suff
 * @param bktSize - bucket size
 * -----------------------------------------*/
void initBkt(BKT *bkt, int bktSize){
    int i;
    for(i = 0; i < bktSize; i++){
        bkt[i].bktList = NULL; 
        bkt[i].bktCnt = 0; 
    }
}


/* -----------------------------------------* 
 * This function can be optimized if higher
 * performance is required.
 *
 * @param kmer - address of of k-mer string
 * @param k - slide window size
 * -----------------------------------------*/
#pragma mta inline
int entryIndex(char *kmer, int k){
    int i;
    int value = 0;

    for(i = 0; i < k; i++){
        value = value*SIGMA + (kmer[i] - 'A');
    }
    
    return value;
}


/* ----------------------------------------------*
 * Bucketing for str <sid>.
 *
 * @param str - seqence itself
 * @param strlen - including ending '$'
 * @param sid - string id. eg. 0, 1, 2, ..., N-1
 * @param bkt - headers for each bucket
 * @param bktSize - bucket size
 * @param sf - mem allocated for all buckets
 * @param k - slide window size
 * ----------------------------------------------*/
#ifdef CRAY_XMT
#pragma mta inline
void slideWindow(char *str, int strLen, int sid, BKT *bkt, int bktSize, SUFFIX *sf, int k){
    int i;
    int stopIndex = strLen - k - 1;
   
    for(i = 0; i <= stopIndex; i++){
        int bktIndex = entryIndex(str+i, k);

        /* prefixed in the bucket list */
        int j = int_fetch_add(&sfIndex, 1);
        sf[j].sid = sid;
        sf[j].pid = i;
        sf[j].next = readfe(&(bkt[bktIndex].bktList));
        writeef(&(bkt[bktIndex].bktList), &sf[j]);
        bkt[bktIndex].bktCnt++;
    }
}

#endif 

#ifndef CRAY_XMT
void slideWindow(char *str, int strLen, int sid, BKT *bkt, int bktSize, SUFFIX *sf, int k){
    int i;
    int stopIndex = strLen - k - 1;
    int bktIndex;
   
    for(i = 0; i <= stopIndex; i++){
        bktIndex = entryIndex(str+i, k);

        /* prefixed in the bucket list */
        sf[sfIndex].sid = sid;
        sf[sfIndex].pid = i;
        sf[sfIndex].next = bkt[bktIndex].bktList;
        bkt[bktIndex].bktList = &sf[sfIndex];

        bkt[bktIndex].bktCnt++;
        sfIndex++;
    }
}
#endif


/* ----------------------------------------------*
 * This function is to bucket all seqs in buckets
 *
 * @param seqs - fasta seqs
 * @param nseqs - #(fasta seqs)
 * @param bkt - global bucket shared by all
 * @param bktSize - bucket size, which is <SIGMA^k>
 * @param suff - 
 * @param suffSize -
 * @param k - slide window size
 * ----------------------------------------------*/
void buildBkt(SEQ *seqs, int nseqs, BKT *bkt, int bktSize, SUFFIX *sf, int sfSize, int k){
#pragma mta assert no alias *seqs, *bkt, *sf
    int i;
    initBkt(bkt, bktSize);
    
    /* slide k-mers for every seqs, and bucket them */
    #pragma mta assert parallel
    for(i = 0; i < nseqs; i++){
        slideWindow(seqs[i].str, seqs[i].strLen, i, bkt, bktSize, sf, k);
    }

    printf("sfIndex=%d\n", sfIndex);
    assert(sfIndex == sfSize);
}


/* ----------------------------------------------*
 * Print bucket bkt[bIndex].
 *
 * @param bkt - bucket pointers
 * @param bIndex - index for each bucket
 * ----------------------------------------------*/
int printBkt(BKT *bkt, int bIndex){
    SUFFIX *p = NULL;
    int i = 0;

    printf("->");
    for(p = bkt[bIndex].bktList; p != NULL; p = p->next){
        printf("[%d, %d]\t", p->sid, p->pid);
        i++;
    }
    printf("\n");
    return i;
}


/* ----------------------------------------------*
 * Print a linked list.
 *
 * @param bkt - bucket pointers
 * @param bIndex - index for each bucket
 * ----------------------------------------------*/
int printBktList(SUFFIX *bktList){
    SUFFIX *p = NULL;
    int i = 0;

    printf("->");
    for(p = bktList; p != NULL; p=p->next){
        printf("[%d, %d]\t", p->sid, p->pid);
        i++;
    }
    printf("\n");
    return i;
}
