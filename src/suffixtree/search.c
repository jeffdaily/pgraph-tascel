#include "search.h"

#pragma mta parallel off
/* -----------------------------------------------------------*
 * sort the bucket lists according to their size. 
 *
 * @param bkt - bkt struct {bktList, bktCnt}
 * @param bktSize - #(bktList)
 * -----------------------------------------------------------*/
void cntSort4Bkt(BKT *bkt, int bktSize){
    int *cnt = NULL;
    int i;
    BKT *nBkt = NULL;
    int maxBktSize = 0;

    /* find the maximum size for bktSize */
    for(i = 0; i < bktSize; i++){
        if(bkt[i].bktCnt > maxBktSize) maxBktSize = bkt[i].bktCnt;
    }
    /* +1 for counters *cnt */
    maxBktSize++;

    cnt = emalloc(maxBktSize*(sizeof *cnt));
    nBkt = emalloc(bktSize*(sizeof *nBkt));

    /* init counter */
    for(i = 0; i < maxBktSize; i++){
        cnt[i] = 0;
    }
    
    /* counting */
    for(i = 0; i < bktSize; i++){
        if(bkt[i].bktCnt >= maxBktSize){
            printf("[Cnt=%d, mSize=%d]\n", bkt[i].bktCnt, maxBktSize); 
            exit(0);
        }
        cnt[bkt[i].bktCnt]++;
    }

    /* prefix sum in a reverse way */
    for(i = maxBktSize-2; i >= 0; i--){
        cnt[i] += cnt[i+1];
    }
    
    /* put the bkt[] into nBkt[] in the descending order */
    for(i = 0; i < bktSize; i++){
        nBkt[cnt[bkt[i].bktCnt]-1] = bkt[i];
        cnt[bkt[i].bktCnt]--;
    }

    /* copy descending order back */
    for(i = 0; i < bktSize; i++){
        bkt[i] = nBkt[i];
    }
    
    free(cnt); 
    free(nBkt); 
}

/* -----------------------------------------------------------*
 * counting sort based on the depth of the stNodes[i].
 * 
 * @param stNodes -
 * @param srtIndex -
 * @param nStNodes -
 * @param maxSeqLen -
 * -----------------------------------------------------------*/
void countSort(STNODE *stNodes, int *srtIndex, int nStNodes, int maxSeqLen){
    int i;
    int *cnt = NULL;
    int depth;

    cnt = emalloc(maxSeqLen*(sizeof *cnt));

    /* init counters */
    for(i = 0; i < maxSeqLen; i++){
        cnt[i] = 0;
    } 

    /* fill in counters */
    for(i = 0; i < nStNodes; i++){
        depth = stNodes[i].depth; /* depth starts from 0 */
        assert(depth < maxSeqLen);
        cnt[depth]++;
    }
    
    /* suffix sum to sort in a reverse way */
    for(i = maxSeqLen-2; i >= 0; i--){
       cnt[i] += cnt[i+1]; 
    }
    
    /* store the sorted index into srtIndex */
    for(i = 0; i < nStNodes; i++){
        depth = stNodes[i].depth;
        srtIndex[cnt[depth]-1] = i;
        cnt[depth]--;
    }

    free(cnt);
}


/* -----------------------------------------------------------*
 * comparison function for qsort().
 * NOTE: it is sorted in a descending way.
 * 
 * @param vp1 -
 * @param vp2 -
 * -----------------------------------------------------------*/
#pragma mta inline
int bktCmp(const void *vp1, const void *vp2){
    BKT *b1 = (BKT *)vp1;
    BKT *b2 = (BKT *)vp2;

    if(b1->bktCnt > b2->bktCnt){
        return -1;
    }else if (b1->bktCnt == b2->bktCnt){
        return 0;
    }else{
        return 1;
    }
}


void printIntArr(int *srtIndex, int size){
    int i;
    for(i = 0; i < size; i++){
        printf("strNdIndex[%ld]=%ld\n", i, srtIndex[i]);
    }
}
