#include "bucket.h"
#include "elib.h"
#include "search.h"
#include "stree.h"

/* -----------------------------------------------------------*
 * sort the bucket lists according to their size. 
 *
 * @param bkt - bkt struct {bktList, size}
 * @param bktSize - #(bktList)
 * -----------------------------------------------------------*/
void cntSort4Bkt(bucket_t *bkt, int bktSize){
    int *cnt = NULL;
    int i;
    bucket_t *nBkt = NULL;
    int maxBktSize = 0;

    /* find the maximum size for bktSize */
    for(i = 0; i < bktSize; i++){
        if(bkt[i].size > maxBktSize) maxBktSize = bkt[i].size;
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
        if(bkt[i].size >= maxBktSize){
            printf("[Cnt=%d, mSize=%d]\n", bkt[i].size, maxBktSize); 
            exit(0);
        }
        cnt[bkt[i].size]++;
    }

    /* prefix sum in a reverse way */
    for(i = maxBktSize-2; i >= 0; i--){
        cnt[i] += cnt[i+1];
    }
    
    /* put the bkt[] into nBkt[] in the descending order */
    for(i = 0; i < bktSize; i++){
        nBkt[cnt[bkt[i].size]-1] = bkt[i];
        cnt[bkt[i].size]--;
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
void countSort(stnode_t *stNodes, int *srtIndex, int nStNodes, int maxSeqLen){
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
int bktCmp(const void *vp1, const void *vp2){
    bucket_t *b1 = (bucket_t *)vp1;
    bucket_t *b2 = (bucket_t *)vp2;

    if(b1->size > b2->size){
        return -1;
    }else if (b1->size == b2->size){
        return 0;
    }else{
        return 1;
    }
}


void printIntArr(int *srtIndex, int size){
    int i;
    for(i = 0; i < size; i++){
        printf("strNdIndex[%d]=%d\n", i, srtIndex[i]);
    }
}
