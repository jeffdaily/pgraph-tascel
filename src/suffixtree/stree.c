/**
 * @file stree.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "bucket.h"
#include "dynamic.h"
#include "elib.h"
#include "stree.h"
#include "timer.h"

int build_tree(
        sequence_t *sequences, size_t n_sequences,
        suffix_t *suffixes, int depth, int window_size,
        stnode_t *st_nodes, int *st_index)
{
    int diffPos;
    int i;
    int j; /* store the index of the internal node */
    suffix_t **heads = NULL;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int hIndex;
    sequence_t *seq = NULL;
    int rLeaf = -1;    

    diffPos = nextDiffPos(sequences, suffixes, depth, window_size);
    if(diffPos == ERROR){ /* can never happen */
        error("wrong when exploring diff chars!");        
        return -1;
    }else if(diffPos == DOL_END){ /* leaf node */
        seq = &sequences[suffixes[0].sid];

        /* depth also includes the '$' */
        st_nodes[*st_index].depth = (seq->strLen - 1) - suffixes[0].pid;
        st_nodes[*st_index].rLeaf = *st_index; /* point to itself */

        compute_lset(suffixes, sequences, st_nodes[*st_index].lset); 

        return (*st_index)++;
    }else{  /* internal node */
        
        j = *st_index; /* store st_index in the stack */
        (*st_index)++; 
        st_nodes[j].depth = diffPos - 1 ;

        heads = st_nodes[j].lset;
        for(i = 0; i < SIGMA; i++) heads[i] = NULL;

        #ifdef DEBUG
        printf("============================\n");
        print_suffixes(suffixes);
        printf("============================\n");
        #endif


        /* partition suffixes into SIGMA sub-buckets */
        for(p = suffixes; p != NULL; p = q){
            q = p->next; 
            hIndex = sequences[p->sid].str[p->pid+diffPos] - 'A';
            assert(hIndex < SIGMA && hIndex >= 0);

            p->next = heads[hIndex];
            heads[hIndex] = p;
        }
        
        #ifdef DEBUG
        for(i = 0; i < SIGMA; i++){
           if(heads[i]) print_suffixes(heads[i]);
        }
        #endif

        /* recursively construct the tree in DFS way */
        for(i = 0; i < SIGMA; i++){
            if(heads[i]){
                /* branching with '$' */
                if(i == (DOLLAR - 'A')){
                    st_nodes[*st_index].depth = (sequences[heads[i]->sid].strLen - 1 - heads[i]->pid);  
                    st_nodes[*st_index].rLeaf = *st_index;

                    compute_lset(heads[i], sequences, st_nodes[*st_index].lset); 
                    rLeaf = (*st_index)++;
                }else{
                    rLeaf = build_tree(sequences, n_sequences, heads[i], diffPos, window_size, st_nodes, st_index);
                }

                /* put it back into NULL */
                heads[i] = NULL; 
            } 
        }

        /* store the right most leaf in the internal node */
        st_nodes[j].rLeaf = rLeaf;
        return st_nodes[j].rLeaf;
    }
}


void compute_lset(suffix_t *suffixes, sequence_t *seqs, suffix_t **lset)
{
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int lIndex;

    for(p = suffixes; p != NULL; p = q){
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
 * @param suffixes -
 * @param blSize -
 * @param window_size -
 * ---------------------------------------------------*/
void procBkt(sequence_t *seqs, int nseqs, suffix_t *suffixes, int blSize, int maxSeqLen, ufind_t *uSet, param_t *param, int ind){

    stnode_t *stNodes = NULL;
    //int *srtIndex = NULL;   /* sorted array based on stnode.depth */
    int depth = param->window_size - 1;
    int stIndex = 0;
    //int *dup = NULL;        /* duplicated entried reduction */

    stNodes = ecalloc(2*blSize, sizeof *stNodes);
    
    build_tree(seqs, nseqs, suffixes, depth, param->window_size, stNodes, &stIndex);

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

void printStnodes(stnode_t *stNodes, int stIndex, int blSize, int ind){

    /* output the tree */
    int i;
    int j;
    FILE *fp = NULL;
    suffix_t *p = NULL;
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
 * @param buckets - bucket
 * @param n_buckets -
 * @param seqs - global fasta seqs
 * @param nseqs -  #seqs in fasta file
 * @param window_size - slide window size
 * @param NN - #chars in all seqs
 * @param nStNodes - #stnodes in constructed forest
 * ---------------------------------------------------*/
void buildForest(bucket_t *buckets, int n_buckets, sequence_t *seqs, int nseqs, int maxSeqLen, ufind_t *uSet, param_t *param){
    int i;
    int sum = 0;
    int cnt = 0;

    for(i = 0; i < n_buckets; i++){
        if(buckets[i].suffixes){
            cnt++;
        }
    }
    printf("cnt=%d\n", cnt);

    for(i = 0; i < n_buckets; i++){
        //double time1 = timer();
        if(buckets[i].suffixes){
            procBkt(seqs, nseqs, buckets[i].suffixes, buckets[i].size, maxSeqLen, uSet, param, i%8192);

            #ifdef DEBUG
            sum += print_suffixes(buckets[i]);
            #endif
        }

        //double time2 = timer();
        //printf("Time for loop <%d> = %lf\n", i, time2 - time1);
    }

    //printf("Time for bucket loop = %lf\n", time2 - time1);

    printf("sum=%d\n", sum);
}

/* ---------------------------------------------------*
 * find the next different position of suffixes.
 * NOTE: if all suffixes end at '$', then report it as
 *       isLeaf = YES
 *
 * return -2: all dollar ending
 *        -1: error, cannot happen
 *         i: next diff position from starting point
 *
 * @param seqs - all fasta seqs
 * @param suffixes - bucket in linked list way 
 * @param depth - depth since root, which has 0 depth
 * @param window_size - slide window size for bucketing
 * ---------------------------------------------------*/
int nextDiffPos(sequence_t *seqs, suffix_t *suffixes, int depth, int window_size){
    int i;
    suffix_t *p = NULL;
    char pCh;
    char cCh;

    assert(suffixes != NULL);
    i = depth; /* position checked until this point */
    p = suffixes;
    
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

        /* reset suffixes to next round of character comparison */
        p = suffixes;

        /* all suffixes ending in '$' */
        if(pCh == DOLLAR) {
            return DOL_END;
        }
    }
}


void printLset(suffix_t **lset){
    int i;

    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    for(i = 0; i < SIGMA; i++){
        if(lset[i]) print_suffixes(lset[i]);
    }
    printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
}


void freeStNodes(stnode_t *stNodes, int size){
    int i;

    for(i = 0L; i < size; i++){
        if(stNodes[i].lset) free(stNodes[i].lset);
    }
    free(stNodes);
}

void freeLset(suffix_t *lset){
    suffix_t *p = NULL;
    suffix_t *q = NULL;

    for(p = lset; p != NULL; p = q){
        q = p->next; 
        free(p);  
    }
}

int lsetSize(suffix_t *lset){
    int size = 0;
    suffix_t *p = NULL;

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

void printStNode(stnode_t *stNodes, int i){
    printf("[%d, %d, %p]\n", stNodes[i].depth, stNodes[i].rLeaf, (void *)stNodes[i].lset);
}

void printStNodes(stnode_t *stNodes, int offset){
    int i;
    //int j;

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
