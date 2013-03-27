#include "pairs.h"
static int generatedPairs = 0;
static int alignedPairs = 0;
static int acceptedPairs = 0;

void printPairs(){
    printf("generatedPairs=%d, alignedPairs=%d, acceptedPairs=%d\n", generatedPairs, alignedPairs, acceptedPairs);
}

/*---------------------------------------------------------------*
 * Generate promising pairs for alignment. 
 *
 * @param stNodes -
 * @param srtIndex -
 * @param nStNodes -
 * @param seqs -
 * @param nSeqs -
 * @param maxSeqLen -
 * @param uSet -
 * @param dup -
 * @param param -
 *---------------------------------------------------------------*/
#pragma mta parallel off
void genPairs(STNODE *stNodes, int *srtIndex, int nStNodes, SEQ *seqs, int nSeqs, int maxSeqLen, UF *uSet, int *dup, PARAM *param){
    int i;
    int j;
    int r;
    STNODE *stnode = NULL;
    int sIndex, eIndex;
    int nps = 0;
    CELL **tbl = NULL;
    int **del = NULL;
    int **ins = NULL;
    int m, n, s, t;
    SUFFIX *p = NULL;
    SUFFIX *q = NULL;
    int f1, f2;
    int s1Len, s2Len;
    int EM;
    int cutOff;         /* cut off value of filter 1 */
    CELL result;

    /* only two rows are allocated */
    assert(NROW == 2);

    EM = param->exactMatchLen;
    cutOff = param->AOL*param->SIM;

    tbl = allocTBL(NROW, maxSeqLen);
    del = allocINT(NROW, maxSeqLen);
    ins = allocINT(NROW, maxSeqLen);


    /* srtIndex maintain an order of NON-increasing depth of stNodes[] */
    for(i = 0; i < nStNodes; i++){
        sIndex = srtIndex[i];
        stnode = &stNodes[sIndex];

        #ifdef DEBUG
        printf("stNode->depth=%d, stnode->rLeaf=%ld, sIndex=%ld\n", stnode->depth, stnode->rLeaf, sIndex);
        #endif

        if(stnode->depth >= EM-1){
            if(stnode->rLeaf == sIndex){ /* leaf node */
                procLeaf(stnode->lset, seqs, nSeqs, tbl, del, ins, uSet, param);
            }else{                       /* internal node */
                eIndex = stnode->rLeaf;

                /* init dup[] for the internal node */
                for(r = 0; r < nSeqs; r++) dup[r] = -1;
                
                /* pairs generation loop for internal node */
                for(m = sIndex+1; m < eIndex; m++){
                    for(n = m+1; n <= eIndex; n++){
                        for(s = 0; s < SIGMA; s++){
                            if(stNodes[m].lset[s]){
                                for(t = 0; t < SIGMA; t++){
                                    if(stNodes[n].lset[t]){
                                       if(s != t){
                                            for(p = stNodes[m].lset[s]; p != NULL; p = p->next){
                                            
                                                /* eliminate pairs */ 
                                                if(dup[p->sid] == -1){
                                                    dup[p->sid] = p->sid;
                                                }else{
                                                    continue;
                                                } 

                                                for(q = stNodes[n].lset[t]; q != NULL; q = q->next){
                                                    f1 = p->sid;
                                                    f2 = q->sid;
                                                    s1Len = seqs[f1].strLen - 1;
                                                    s2Len = seqs[f2].strLen - 1;
                                                    
                                                    if(f1 == f2) continue;
                                                    continue;

                                                    generatedPairs++;

                                                    if(find(uSet, f1) != find(uSet, f2)){
                                                        alignedPairs++;
                                                        if(s1Len <= s2Len){
                                                            if(100*s1Len < cutOff*s2Len) continue;
                                                            affineGapAlign(seqs[f1].str, s1Len,
                                                                           seqs[f2].str, s2Len,
                                                                           &result, tbl, del, ins);
                                                        }else{
                                                            if(100*s2Len < cutOff*s1Len) continue;
                                                            affineGapAlign(seqs[f2].str, s2Len,
                                                                           seqs[f1].str, s1Len,
                                                                           &result, tbl, del, ins);

                                                        }

                                                        /* check if it is an edge */
                                                        if(isEdge(&result, seqs[f1].str, s1Len, seqs[f2].str, s2Len, param)){
                                                            acceptedPairs++;
                                                            union_elems(uSet, f1, f2);
                                                        }
                                                        //printf("PAIR:[%d, %d]\n", f1, f2);
                                                    }
                                                }
                                            }
                                        } 
                                    }
                                }
                            }
                        }
                    }
                }
    



                /* merge the lsets of subtree */ 
                for(m = 0; m < SIGMA; m++){
                    p = NULL;
                    for(j = sIndex+1; j <= eIndex; j++){
                        if((q = stNodes[j].lset[m])){

                            /* empty the subtree's ptrs array */
                            stNodes[j].lset[m] = NULL;
                            if(p == NULL) {
                                p = q;
                                stNodes[sIndex].lset[m] = q;
                            } else p->next = q;

                            /* walk to the end */
                            while(p->next) p = p->next;
                        }  
                    }
                }
            }
        }else{
            /* stnodes are sorted, so later part 
             * will not satisfy EM cutoff */        
            break; 
        }
    }  

    /* free */
    freeTBL(tbl, NROW);
    freeINT(del, NROW);
    freeINT(ins, NROW);
}

/*---------------------------------------------------------------*
 * This function implements the pair generation algorithm for leaf
 * nodes. 
 *
 *    BEGIN - intra/inter cross. O/W - intra cross. 
 *
 * @param lset -
 * @param seqs -
 * @param nSeqs -
 * @param tbl -
 * @param ins -
 * @param del -
 *---------------------------------------------------------------*/
#pragma mta parallel off
void procLeaf(SUFFIX **lset, SEQ *seqs, int nSeqs, CELL **tbl, int **ins, int **del, UF *uSet, PARAM *param){
    int i;
    int j;
    SUFFIX *p = NULL;
    SUFFIX *q = NULL;
    int f1, f2;
    int s1Len, s2Len;
    CELL result;
    int cutOff;

    cutOff = param->AOL*param->SIM;

    for(i = 0; i < SIGMA; i++){
        if(lset[i]){
            if(i == BEGIN - 'A'){ /* inter cross */
                for(p = lset[i]; p != NULL; p = p->next){
                    for(q = p->next; q != NULL; q = q->next){
                        f1 = p->sid;
                        f2 = q->sid;

                        if(f1 == f2) continue;
                        continue;

                        generatedPairs++;

                        s1Len = seqs[f1].strLen - 1;
                        s2Len = seqs[f2].strLen - 1;

                        if(find(uSet, f1) != find(uSet, f2)){
                            alignedPairs++;
                            if(s1Len <= s2Len){
                                if(100*s1Len < cutOff*s2Len) continue;
                                affineGapAlign(seqs[f1].str, s1Len,
                                               seqs[f2].str, s2Len,
                                               &result, tbl, del, ins);
                            }else{
                                if(100*s2Len < cutOff*s1Len) continue;
                                affineGapAlign(seqs[f2].str, s2Len,
                                               seqs[f1].str, s1Len,
                                               &result, tbl, del, ins);

                            }

                            /* check if it is an edge */
                           if(isEdge(&result, seqs[f1].str, s1Len, seqs[f2].str, s2Len, param)){
                                acceptedPairs++;
                                union_elems(uSet, f1, f2);
                           }
                           //printf("PAIR:[%d, %d]\n", f1, f2);
                        }
                        //printf("[%d, %d] - score: %d\n", f1, f2, result.score);
                    }
                } 
            }        
           
            /* intra cross */
            for(j = i+1; j < SIGMA; j++){
                if(lset[j]){
                    for(p = lset[i]; p != NULL; p = p->next){
                        for(q = lset[j]; q != NULL; q = q->next){
                            f1 = p->sid;
                            f2 = q->sid;

                            if(f1 == f2) continue;
                            continue;
                            generatedPairs++;

                            s1Len = seqs[f1].strLen - 1;
                            s2Len = seqs[f2].strLen - 1;

                            if(find(uSet, f1) != find(uSet, f2)){
                                alignedPairs++;
                                if(s1Len <= s2Len){
                                    if(100*s1Len < cutOff*s2Len) continue;
                                    affineGapAlign(seqs[f1].str, s1Len,
                                                   seqs[f2].str, s2Len,
                                                   &result, tbl, del, ins);
                                }else{
                                    if(100*s2Len < cutOff*s1Len) continue;
                                    affineGapAlign(seqs[f2].str, s2Len,
                                                   seqs[f1].str, s1Len,
                                                   &result, tbl, del, ins);

                                }

                                /* check if it is an edge */
                                if(isEdge(&result, seqs[f1].str, s1Len, seqs[f2].str, s2Len, param)){
                                    acceptedPairs++;
                                    union_elems(uSet, f1, f2);
                                }
                                //printf("PAIR:[%d, %d]\n", f1, f2);
                            }
                            //printf("[%d, %d] - score: %d\n", f1, f2, result.score);
                        }
                    }
                }
            }
        }
    }    
}



#pragma mta inline
int isEdge(CELL *result, char *s1, int s1Len, char *s2, int s2Len, PARAM *param){
    int sscore;
    int maxLen;
    int nmatch;
    
    if(result->score <= 0) return FALSE;
    
    /* DO NOT need to compute the sscore value every time, if the later check
     * failed at the first step, the sscore computation is wasted */
    if(s1Len > s2Len){
        maxLen = s1Len;
        sscore = selfScore(s1, s1Len);
    }else{
        maxLen = s2Len;
        sscore = selfScore(s2, s2Len);
    }   
    
    nmatch = result->ndig;
    
    /* order the condition in strict->loose way, performance perspective 
     * comparison using integers, no overflow could happen */
    if((10*result->alen >= param->AOL*maxLen) 
        && (10*nmatch >= param->SIM*result->alen)
        && (10*result->score >= param->OS*sscore)){
        return TRUE;
    }else{
        return FALSE;
    }   
}

