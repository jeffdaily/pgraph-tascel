#include "stree.h"
#include "timer.h"

/* -----------------------------------------------------------*
 * build tree for each bucket, each bucket will be
 * branched into SIGMA different sub-buckets in a DFS
 * way, the total space will be (26*4*depth) bytes
 *
 * NOTE: bktList is modified after this function, which 
 * also means the previous bucket does not exist at all
 *
 * @param seqs - all fasta seqs
 * @param bktList - suffixes list, modified thereafter
 * @param depth - characters compared since root, which has
 *                0 depth
 * @param k - slide window size
 * @param stNodes - stnodes array in DFS order
 * -----------------------------------------------------------*/
#pragma mta parallel off
int buildTree(SEQ *seqs, int nseqs, SUFFIX *bktList, int depth, int k, STNODE *stNodes, int *stIndex){
    int diffPos;
    int i;
    int j; /* store the index of the internal node */
    SUFFIX **heads = NULL;
    SUFFIX *p = NULL;
    SUFFIX *q = NULL;
    int hIndex;
    SEQ *seq = NULL;
    int rLeaf = -1;    

    diffPos = nextDiffPos(seqs, bktList, depth, k);
    if(diffPos == ERROR){ /* can never happen */
        error("wrong when exploring diff chars!");        
        return -1;
    }else if(diffPos == DOL_END){ /* leaf node */
        seq = &seqs[bktList[0].sid];

        /* depth also includes the '$' */
        stNodes[*stIndex].depth = (seq->strLen - 1) - bktList[0].pid;
        stNodes[*stIndex].rLeaf = *stIndex; /* point to itself */

        compLset(bktList, seqs, nseqs, stNodes[*stIndex].lset); 

        return (*stIndex)++;
    }else{  /* internal node */
        
        j = *stIndex; /* store stIndex in the stack */
        (*stIndex)++; 
        stNodes[j].depth = diffPos - 1 ;

        heads = stNodes[j].lset;
        for(i = 0; i < SIGMA; i++) heads[i] = NULL;

        #ifdef DEBUG
        printf("============================\n");
        printBktList(bktList);
        printf("============================\n");
        #endif


        /* partition bktList into SIGMA sub-buckets */
        for(p = bktList; p != NULL; p = q){
            q = p->next; 
            hIndex = seqs[p->sid].str[p->pid+diffPos] - 'A';
            assert(hIndex < SIGMA && hIndex >= 0);

            p->next = heads[hIndex];
            heads[hIndex] = p;
        }
        
        #ifdef DEBUG
        for(i = 0; i < SIGMA; i++){
           if(heads[i]) printBktList(heads[i]);
        }
        #endif

        /* recursively construct the tree in DFS way */
        for(i = 0; i < SIGMA; i++){
            if(heads[i]){
                /* branching with '$' */
                if(i == (DOLLAR - 'A')){
                    stNodes[*stIndex].depth = (seqs[heads[i]->sid].strLen - 1 - heads[i]->pid);  
                    stNodes[*stIndex].rLeaf = *stIndex;

                    compLset(heads[i], seqs, nseqs, stNodes[*stIndex].lset); 
                    rLeaf = (*stIndex)++;
                }else{
                    rLeaf = buildTree(seqs, nseqs, heads[i], diffPos, k, stNodes, stIndex);
                }

                /* put it back into NULL */
                heads[i] = NULL; 
            } 
        }

        /* store the right most leaf in the internal node */
        stNodes[j].rLeaf = rLeaf;
        return stNodes[j].rLeaf;
    }
}


/* --------------------------------------------------------*
 * compute lset according to the left character of pid to 
 * ensure the left maximality. Mark the special case of
 * whole sequence. Use BEGIN='O' as the special character.
 *
 *
 * @param bktList - 
 * @param seqs - all seqs info
 * @param nseqs - #seqs
 * @param lset - 
 * --------------------------------------------------------*/
void compLset(SUFFIX *bktList, SEQ *seqs, int nseqs, SUFFIX **lset){
    SUFFIX *p = NULL;
    SUFFIX *q = NULL;
    int lIndex;

    for(p = bktList; p != NULL; p = q){
        q = p->next;
        if(p->pid == 0){
            lIndex = BEGIN - 'A';
            p->next = lset[lIndex];
            lset[lIndex] = p; 
        }else{
            lIndex = seqs[p->sid].str[p->pid-1] - 'A';
            p->next = lset[lIndex];
            lset[lIndex] = p;
        }
    }
}


/* ---------------------------------------------------*
 * This function constructs tree, sort it, and generate
 * pairs.
 * 
 * @param seqs -
 * @param nseqs -
 * @param bktList -
 * @param blSize -
 * @param k -
 * ---------------------------------------------------*/
void procBkt(SEQ *seqs, int nseqs, SUFFIX *bktList, int blSize, int maxSeqLen, UF *uSet, PARAM *param, int ind){

    STNODE *stNodes = NULL;
    int *srtIndex = NULL;   /* sorted array based on stnode.depth */
    int depth = param->k - 1;
    int stIndex = 0;
    int *dup = NULL;        /* duplicated entried reduction */

    stNodes = ecalloc(2*blSize, sizeof *stNodes);
    
    buildTree(seqs, nseqs, bktList, depth, param->k, stNodes, &stIndex);

    /* output stnodes into file */
    printStnodes(stNodes, stIndex, blSize, ind);

    /* end of outputing tree */
    /* sort the stnodes in another array */
    //srtIndex = emalloc(stIndex*(sizeof *srtIndex));
    //countSort(stNodes, srtIndex, stIndex, maxSeqLen);
    

    /* pairs generation and alignment */
    //dup = emalloc(nseqs*(sizeof *dup));
    //genPairs(stNodes, srtIndex, stIndex, seqs, nseqs, maxSeqLen, uSet, dup, param);  
    

    //free(dup);
    //free(srtIndex);
    free(stNodes);
}

void printStnodes(STNODE *stNodes, int stIndex, int blSize, int ind){

    /* output the tree */
    int i;
    int j;
    FILE *fp = NULL;
    SUFFIX *p = NULL;
    int m;
    char frFile[200];


    sprintf(frFile, "./320k_tree/forest_%d", ind);
    fp = efopen(frFile, "a+");
    fprintf(fp, "<stSize: %d, bktSize: %d\n", stIndex, blSize);
    for(i = 0; i< stIndex; i++){
        fprintf(fp, "=%d, %d", stNodes[i].depth, stNodes[i].rLeaf);    
        for(m = 0, j = 0; j < SIGMA; j++){
            //fprintf(fp, "%d", j);
            if(stNodes[i].lset[j]){
                fprintf(fp, "\n");
                m++;
                fprintf(fp, "%d ", j);
                for(p = stNodes[i].lset[j]; p != NULL; p = p->next){
                    fprintf(fp, "[%d,%d] ", p->sid, p->pid);
                }
            }
        }
        fprintf(fp, "\n");

    }
    fclose(fp);
}


/* ---------------------------------------------------*
 * build a tree for each bucket
 *
 * @param bkt - bucket
 * @param bktSize -
 * @param seqs - global fasta seqs
 * @param nseqs -  #seqs in fasta file
 * @param k - slide window size
 * @param NN - #chars in all seqs
 * @param nStNodes - #stnodes in constructed forest
 * ---------------------------------------------------*/
#pragma mta parallel on
void buildForest(BKT *bkt, int bktSize, SEQ *seqs, int nseqs, int maxSeqLen, UF *uSet, PARAM *param){
    int i;
    int sum = 0;
    int ind = 0;
    int cnt = 0;

    for(i = 0; i < bktSize; i++){
        if(bkt[i].bktList){
            cnt++;
        }
    }
    printf("cnt=%d\n", cnt);

    #pragma mta assert parallel
    #pragma mta use 100 streams
    #pragma mta dynamic schedule
    for(i = 0; i < bktSize; i++){
        //double time1 = timer();
        if(bkt[i].bktList){
            procBkt(seqs, nseqs, bkt[i].bktList, bkt[i].bktCnt, maxSeqLen, uSet, param, i%8192);

            #ifdef DEBUG
            sum += printBktList(bkt[i]);
            #endif
        }

        //double time2 = timer();
        //printf("Time for loop <%d> = %lf\n", i, time2 - time1);
    }

    //printf("Time for bucket loop = %lf\n", time2 - time1);

    printf("sum=%d\n", sum);
}

/* ---------------------------------------------------*
 * find the next different position of bktList.
 * NOTE: if all suffixes end at '$', then report it as
 *       isLeaf = YES
 *
 * return -2: all dollar ending
 *        -1: error, cannot happen
 *         i: next diff position from starting point
 *
 * @param seqs - all fasta seqs
 * @param bktList - bucket in linked list way 
 * @param depth - depth since root, which has 0 depth
 * @param k - slide window size for bucketing
 * ---------------------------------------------------*/
#pragma mta parallel off 
int nextDiffPos(SEQ *seqs, SUFFIX *bktList, int depth, int k){
    int i;
    SUFFIX *p = NULL;
    char pCh;
    char cCh;

    assert(bktList != NULL);
    i = depth; /* position checked until this point */
    p = bktList;
    
    while(1){
        i++; /* step forward one more char for comparison */
        assert((i+p->pid) < seqs[p->sid].strLen);
        pCh = seqs[p->sid].str[p->pid+i];

        for(p = p->next; p != NULL; p = p->next){
            assert((i+p->pid) < seqs[p->sid].strLen);
            cCh = seqs[p->sid].str[p->pid+i];
            if(cCh != pCh){
                return i;
            }
        }

        /* reset bktList to next round of character comparison */
        p = bktList;

        /* all suffixes ending in '$' */
        if(pCh == DOLLAR) {
            return DOL_END;
        }
    }
}


void printLset(SUFFIX **lset){
    int i;

    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    for(i = 0; i < SIGMA; i++){
        if(lset[i]) printBktList(lset[i]);
    }
    printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
}


void freeStNodes(STNODE *stNodes, int size){
    int i;

    for(i = 0L; i < size; i++){
        if(stNodes[i].lset) free(stNodes[i].lset);
    }
    free(stNodes);
}

void freeLset(SUFFIX *lset){
    SUFFIX *p = NULL;
    SUFFIX *q = NULL;

    for(p = lset; p != NULL; p = q){
        q = p->next; 
        free(p);  
    }
}

int lsetSize(SUFFIX *lset){
    int size = 0;
    SUFFIX *p = NULL;

    for(p = lset; p != NULL; p = p->next){
        size++;
    }
    return size;
}

void printCnt(int *cnt, int size){
    int i;
    
    for(i = 0; i < size; i++){
        printf("cnt[%d]=%d\n", i, cnt[i]); 
    }
}

void printStNode(STNODE *stNodes, int i){
    printf("[%d, %d, %p]\n", stNodes[i].depth, stNodes[i].rLeaf, (void *)stNodes[i].lset);
}

void printStNodes(STNODE *stNodes, int offset){
    int i;
    int j;

    for(i = 0; i < offset; i++){
        printf("[%d, %d]\n", stNodes[i].depth, stNodes[i].rLeaf);
        printf("-------------------------------\n");
        /*for(j = 0; j < SIGMA; j++){
            printf("%p\n", (void *)stNodes[i].lset[j]);
        }*/
        //printf("-------------------------------\n");
    }
    printf("\n");
}
