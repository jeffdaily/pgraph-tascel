#include "pairs.h"
static int generatedPairs = 0;
static int alignedPairs = 0;
static int acceptedPairs = 0;

void print_pairs()
{
    printf("generatedPairs=%d, alignedPairs=%d, acceptedPairs=%d\n",
            generatedPairs, alignedPairs, acceptedPairs);
}


void genPairs(stnode_t *stNodes, int *srtIndex, int nStNodes, sequence_t *seqs, int nSeqs, int maxSeqLen, ufind_t *uSet, int *dup, param_t *param)
{
    int i;
    int j;
    int r;
    stnode_t *stnode = NULL;
    int sIndex, eIndex;
    cell_t **tbl = NULL;
    int **del = NULL;
    int **ins = NULL;
    int n;
    size_t m, s, t;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int f1, f2;
    int s1Len, s2Len;
    int EM;
    int cutOff;         /* cut off value of filter 1 */
    cell_t result;

    /* only two rows are allocated */
    assert(NROW == 2);

    EM = param->exact_match_len;
    cutOff = param->AOL*param->SIM;

    tbl = alloc_tbl(NROW, maxSeqLen);
    del = alloc_int(NROW, maxSeqLen);
    ins = alloc_int(NROW, maxSeqLen);


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
                                                            affine_gap_align(seqs[f1].str, s1Len,
                                                                    seqs[f2].str, s2Len,
                                                                    &result, tbl, del, ins);
                                                        }else{
                                                            if(100*s2Len < cutOff*s1Len) continue;
                                                            affine_gap_align(seqs[f2].str, s2Len,
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
    free_tbl(tbl, NROW);
    free_int(del, NROW);
    free_int(ins, NROW);
}


void procLeaf(suffix_t **lset, sequence_t *seqs, int nSeqs, cell_t **tbl, int **ins, int **del, ufind_t *uSet, param_t *param)
{
    size_t i;
    size_t j;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int f1, f2;
    int s1Len, s2Len;
    cell_t result;
    int cutOff;

    cutOff = param->AOL*param->SIM;

    for(i = 0; i < SIGMA; i++){
        if(lset[i]){
            if(i == BEGIN - 'A'){ /* inter cross */
                for(p = lset[i]; p != NULL; p = p->next){
                    for(q = p->next; q != NULL; q = q->next){
                        f1 = p->sid;
                        f2 = q->sid;
                        assert(f1 < nSeqs);
                        assert(f2 < nSeqs);

                        if(f1 == f2) continue;
                        continue;

                        generatedPairs++;

                        s1Len = seqs[f1].strLen - 1;
                        s2Len = seqs[f2].strLen - 1;

                        if(find(uSet, f1) != find(uSet, f2)){
                            alignedPairs++;
                            if(s1Len <= s2Len){
                                if(100*s1Len < cutOff*s2Len) continue;
                                affine_gap_align(seqs[f1].str, s1Len,
                                               seqs[f2].str, s2Len,
                                               &result, tbl, del, ins);
                            }else{
                                if(100*s2Len < cutOff*s1Len) continue;
                                affine_gap_align(seqs[f2].str, s2Len,
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
                                    affine_gap_align(seqs[f1].str, s1Len,
                                                   seqs[f2].str, s2Len,
                                                   &result, tbl, del, ins);
                                }else{
                                    if(100*s2Len < cutOff*s1Len) continue;
                                    affine_gap_align(seqs[f2].str, s2Len,
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



int isEdge(cell_t *result, char *s1, int s1Len, char *s2, int s2Len, param_t *param){
    int sscore;
    int maxLen;
    int nmatch;
    
    if(result->score <= 0) return FALSE;
    
    /* DO NOT need to compute the sscore value every time, if the later check
     * failed at the first step, the sscore computation is wasted */
    if(s1Len > s2Len){
        maxLen = s1Len;
        sscore = self_score(s1, s1Len);
    }else{
        maxLen = s2Len;
        sscore = self_score(s2, s2Len);
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

