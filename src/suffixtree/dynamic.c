#include "dynamic.h"

#pragma mta parallel off
static int blosum62[24][24] = { 
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4}, 
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4}, 
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4}, 
{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4}, 
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4}, 
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4}, 
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4}, 
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4}, 
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4}, 
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4}, 
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4}, 
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4}, 
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4}, 
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4}, 
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4}, 
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4}, 
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4}, 
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4}, 
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4}, 
{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4}, 
{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4}, 
{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4}, 
{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};


/* J stands for '*', 'O' and 'U' are dummy ones to make 26 */
static char AA[]={'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
                  'B', 'Z', 'X', 'J', 'O', 'U'};

static int map[SIGMA];

void initMap(int nAA){
    int i;
    for(i = 0; i < nAA; i++){
        map[AA[i] - 'A'] = i;
    }
}

#pragma mta inline
int selfScore(char *s, int ns){
    int i;
    int j;
    int score = 0;

    for(i = 0; i < ns; i++){
        j = map[s[i] - 'A']; 
        score += blosum62[j][j];
    }

    return score;
}

/*--------------------------------------------*
 * alloc 2-dimension CELL array.
 * @param nrow - 
 * @param ncol - strlen(s2), alloc ONE more
 *               dor dynamic align.
 *--------------------------------------------*/
CELL **allocTBL(int nrow, int ncol){
    int i;
    CELL **tbl = NULL;

    tbl = emalloc(nrow*(sizeof *tbl));
    for(i = 0; i < nrow; i++){
        /* +1 for dynamic align */
        tbl[i] = emalloc((ncol+1)*(sizeof *tbl[i]));
    }

    return tbl;
}


/*--------------------------------------------*
 * alloc 2-dimension int array.
 * @param nrow - 
 * @param ncol - strlen(s2), alloc ONE more
 *               dor dynamic align.
 *--------------------------------------------*/
int **allocINT(int nrow, int ncol){
    int **tbl = NULL;
    int i;
    
    tbl = emalloc(nrow*(sizeof *tbl));
    for(i = 0; i < nrow; i++){
        /* +1 for dynamic align */
        tbl[i] = emalloc((ncol+1)*(sizeof *tbl[i]));
    }

    return tbl;
}


/* free function, symetric to allocTBL() */
void freeTBL(CELL **tbl, int nrow){
    int i;
    for(i = 0; i < nrow; i++){
        free(tbl[i]);
    }
    free(tbl);
}


/* free function, symetric to allocINT() */
void freeINT(int **tbl, int nrow){
    int i;
    for(i = 0; i < nrow; i++){
        free(tbl[i]);
    }
    free(tbl);
}


/*------------------------------------------------------------*
 * Implementation of affine gap pairwise sequence alignment, it
 * is a space efficient version: only two rows are required; also
 * mem for all dynamic tables are allocated ONLY ONCE.
 *
 * @param s1 - sequence s1
 * @param s1Len - sequence length of <s1>, strlen(s1)
 * @param s2 - sequence s2
 * @param s2Len - sequence length of <s2>, strlen(s2)
 * @param result - alignment result <score, ndig, alen>
 * @param tbl - pre-allocated score table
 * @param del - pre-allocated deletion table
 * @param ins - pre-allocated insertion table
 *------------------------------------------------------------*/
#pragma mta inline
void ffineGapAlign(char *s1, int s1Len, char *s2, int s2Len, CELL *result, CELL **tbl, int **del, int **ins){
    int i, j;
    int cr, pr;
    int dig, up, left, maxScore;
    char ch1, ch2;  /* character s[i] and s[j] */

    /* struct can be ONLY initialized as it is declared ??*/
    CELL lastCol = {INT_MIN, 0, 0};
    CELL lastRow = {INT_MIN, 0, 0};

    //assert(s1Len>0 && s2Len>0);

    cr = 1;
    pr = 0;

    /* init first row of 3 tables */
    tbl[0][0].score = 0;
    tbl[0][0].ndig = 0;
    tbl[0][0].alen = 0;
    del[0][0] = 0;
    ins[0][0] = 0;

    for(j = 1; j <= s2Len; j++){
        tbl[0][j].score = OPEN + j*GAP;
        tbl[0][j].ndig = 0;
        tbl[0][j].alen = 0;

        del[0][j] = INT_MIN;
        ins[0][j] = OPEN + j*GAP;
    }
   
    for(i = 1; i <= s1Len; i++){
        ch1 = s1[i-1];
        cr = CROW(i);
        pr = PROW(i); 

        /* init first column of 3 tables */
        tbl[cr][0].score = OPEN + i*GAP;
        tbl[cr][0].ndig = 0;
        tbl[cr][0].alen = 0;

        del[cr][0] = OPEN + i*GAP;
        ins[cr][0] = INT_MIN;

        for(j = 1; j <= s2Len; j++){
            ch2 = s2[j-1];

            /* overflow could happen, INT_MIN-1 = 2147483647
             * #define NEG_ADD(x, y) \
             *     (((y)<0)&&((x)<(INT_MIN-y)) ? INT_MIN : (x)+(y)) */
            up = MAX(tbl[pr][j].score+OPEN+GAP, NEG_ADD(del[pr][j], GAP)); 
            del[cr][j] = up;
            left = MAX(tbl[cr][j-1].score+OPEN+GAP, NEG_ADD(ins[cr][j-1], GAP));  
            ins[cr][j] = left;
            maxScore = (up >= left)? up : left;
            
            /* blosum62[map[ch1-'A']][map[ch2-'A']]; */
            dig = tbl[pr][j-1].score + BLOSUM62(map, ch1, ch2); 
            if(dig >= maxScore) maxScore = dig;
            tbl[cr][j].score = maxScore;

            #ifdef DEBUG
            printf("up=%d, left=%d, dig=%d, <%c,%c>\n", up, left, dig, ch1, ch2);
            #endif

            if(maxScore == dig){
                tbl[cr][j].ndig = tbl[pr][j-1].ndig + ((ch1 == ch2) ? 1 : 0);
                tbl[cr][j].alen = tbl[pr][j-1].alen + 1;
            }else if (maxScore == up){
                tbl[cr][j].ndig = tbl[pr][j].ndig;
                tbl[cr][j].alen = tbl[pr][j].alen + 1;
            }else{
                tbl[cr][j].ndig = tbl[cr][j-1].ndig;
                tbl[cr][j].alen = tbl[cr][j-1].alen + 1;
            }
            
            /* track of the maximum last row */
            if(i == s1Len){
                lastRow = (tbl[cr][j].score > lastRow.score) ? tbl[cr][j] : lastRow;
            }
        }

        //assert(j == (s2Len+1));
        /* update the maximum of last column */
        lastCol = (tbl[cr][s2Len].score > lastCol.score)? tbl[cr][s2Len] : lastCol; 
    }

    *result = (lastCol.score > lastRow.score) ? lastCol : lastRow;
}


/*------------------------------------------------------------*
 * Implementation of affine gap pairwise sequence alignment, it
 * is a space efficient version: only two rows are required; also
 * mem for all dynamic tables are allocated ONLY ONCE.
 *
 * @param s1 - sequence s1
 * @param s1Len - sequence length of <s1>, strlen(s1)
 * @param s2 - sequence s2
 * @param s2Len - sequence length of <s2>, strlen(s2)
 * @param result - alignment result <score, ndig, alen>
 * @param tbl - pre-allocated score table
 * @param del - pre-allocated deletion table
 * @param ins - pre-allocated insertion table
 *------------------------------------------------------------*/
void affineGapAlign(char *s1, int s1Len, char *s2, int s2Len, CELL *result, CELL **tbl, int **del, int **ins){
    int i, j;
    int maxScore;
    CELL * maxRecord;

    CELL * tblCR;
    CELL * tblPR;
    int  * delCR;
    int  * delPR;

    tblCR = tbl[0];
    delCR = del[0];

    /* init first row of 3 tables */
    tblCR[0].score = 0;
    tblCR[0].ndig = 0;
    tblCR[0].alen = 0;
    delCR[0] = 0;
    
    for (j = 1; j <= s2Len; j++){
        tblCR[j].score = OPEN + j * GAP;
        tblCR[j].ndig = 0;
        tblCR[j].alen = 0;

        delCR[j] = MIN_VAL;
    }

    maxScore  = tblCR[s2Len].score;
    maxRecord = tblCR + s2Len;

    for (i = 1; i <= s1Len; i++){
        char ch1 = s1[i-1];
        int * BlosumRow = blosum62[map[ch1 - 'A']];

        int cr = CROW(i);
        int pr = PROW(i);
        tblCR  = tbl[cr];
        tblPR  = tbl[pr];
        delCR  = del[cr];
        delPR  = del[pr];

        /* j = 0 */
        int ins    = MIN_VAL;

        int Nscore = tblPR[0].score;
        int Nndig  = tblPR[0].ndig;
        int Nalen  = tblPR[0].alen;

        int Wscore = OPEN + i * GAP;
        int Wndig  = 0;
        int Walen  = 0;

        delCR[0]       = OPEN + i * GAP;
        tblCR[0].score = Wscore;
        tblCR[0].ndig  = Wndig;
        tblCR[0].alen  = Walen;

        for (j = 1; j <= s2Len; j++){
            char ch2    = s2[j-1];
            int NWscore = Nscore;
            int NWndig  = Nndig;
            int NWalen  = Nalen;

            Nscore = tblPR[j].score;
            Nndig  = tblPR[j].ndig;
            Nalen  = tblPR[j].alen;

            int up   = MAX(Nscore + OPEN, delPR[j]) + GAP;
            int left = MAX(Wscore + OPEN, ins     ) + GAP;
            int dig  = NWscore + BlosumRow[map[ch2 - 'A']];

            delCR[j] = up;
            ins      = left;
            
            if ((dig >= up) && (dig >= left)) {
                Wscore = dig;
                Wndig  = NWndig + (ch1 == ch2);
                Walen  = NWalen + 1;
            } else if (up >= left) {
                Wscore = up;
                Wndig  = Nndig;
                Walen  = Nalen + 1;
            } else {
                Wscore = left;
                Wndig  = Wndig;
                Walen  = Walen + 1;
            }

            tblCR[j].score = Wscore;
            tblCR[j].ndig  = Wndig;
            tblCR[j].alen  = Walen;
        }

        if (Wscore > maxScore) {maxScore = Wscore; maxRecord = tblCR + s2Len;}
 
    } /* end of i loop */

    tblCR = tbl[ CROW(s1Len) ];
    for (j = 0; j <= s2Len; j++)
        if (tblCR[j].score > maxScore) {maxScore = tblCR[j].score; maxRecord = tblCR + j;}
            
    * result = * maxRecord;
}


void printRow(CELL **tbl, int i, int ncol){
    int j;
    printf("[");
    for(j = 0; j <= ncol; j++){
        printf("%d, %d\t", tbl[i][j].score, tbl[i][j].ndig);
    }
    printf("]\n");
}
