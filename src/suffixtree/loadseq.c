#include "loadseq.h"

#ifdef CRAY_XMT

#include <sys/stat.h>
#include <luc/luc_exported.h>
#include <snapshot/client.h>

/*------------------------------------------------*
 * This function reads fasta file into SEQ *seqs.
 * NOTE: a '$' signed is appended at the end 
 *
 * @param fileName - fasta file name
 * @param nSeqs - number of seqs in fasta file
 * @param seqs - all fasta seqs
 * @param nChars - #AAs in fasta seqs, including '$'
 * @param maxSeqLen - maximal seq length
 *------------------------------------------------*/
void loadAllSeqs(char *fileName, int nSeqs, SEQ *seqs, int *nChars, int *maxSeqLen){
    char *buffer = NULL;
    struct stat buf;
    luc_error_t err;
    int snap_err = 0;
    int maxLen = 0;

    stat(fileName, &buf);
    int buffer_size = buf.st_size;

    buffer = (char *) malloc((buffer_size + 1) * sizeof(char));
    if (buffer == NULL) {printf("File buffer not allocated\n"); exit(1);}

    err = snap_init();
    if (err != SNAP_ERR_OK) {printf("snap_init failed\n"); exit(1);}

    err = snap_restore(fileName, buffer, buffer_size, &snap_err);
    if (err != SNAP_ERR_OK) {printf("snap_restore failed\n"); exit(1);}

    int NN = 0; 
    int next_flag = 0;
    int *flags = (int *)malloc(nSeqs * sizeof(int));

    /* turn the parallel off to keep the order of sequences */
    #pragma mta assert no dependence
    for (int i = 0; i < buffer_size; i++)
        if (buffer[i] == FASTA_FLAG) flags[ int_fetch_add(&next_flag, 1) ] = i;


    #pragma mta assert no dependence
    for (int i = 0; i < nSeqs; i++) {
        
        int ndx = flags[i] + 1;
        seqs[i].gid = buffer + ndx;

        while (!isspace(buffer[ndx])) ndx++;
        buffer[ndx] = '\0';

        ndx++;
        int strLen  = 1;
        seqs[i].str = buffer + ndx;
        while (!isspace(buffer[ndx])) {ndx++; strLen++;}

        buffer[ndx]     = DOLLAR;
        buffer[ndx + 1] = '\0';
        seqs[i].strLen  = strLen;
        if(strLen > maxLen) maxLen = strLen;
        NN             += strLen;
    }
    
    *nChars = NN;
    *maxSeqLen = maxLen;    
    
    free(flags);
}

#endif

#ifndef CRAY_XMT
void loadAllSeqs(char *fileName, int nSeqs, SEQ *seqs, int *nChars, int *maxSeqLen){
    FILE *fp = NULL;
    char line[MAX_LINE_LEN];
    int lineLen;
    int i = 0; 

    fp = efopen(fileName, "r"); 

    *nChars = 0L;
    *maxSeqLen = 0;
    while(fgets(line, MAX_LINE_LEN, fp)){
        lineLen = strlen(line); 
        assert(lineLen <= MAX_LINE_LEN);
        /* add an '$' end to each seq */
        line[lineLen-1] = '\0';

        if(line[0] == FASTA_FLAG){
            seqs[i].gid = estrdup(line);
        }else if(isalpha(line[0])){
            line[lineLen-1] = DOLLAR;
            line[lineLen] = '\0';

            /* strlen does not include the '\0' */
            seqs[i].strLen = lineLen;
            *nChars += lineLen;
            if(lineLen > *maxSeqLen){
                *maxSeqLen = lineLen;
            }

            seqs[i++].str = estrdup(line); 
        }else{
            warn("empty line in fasta file? will continue!!"); 
        }
    }

    fclose(fp);
    assert(i == nSeqs);
}
#endif


void printSeq(SEQ *seqs, int index){
    printf("StrLen=%d\n", seqs[index].strLen);
    printf("--------------------\n");
    printf("%s\n%s\n", seqs[index].gid, seqs[index].str);
}


int freeSeqs(SEQ *seqs, int nSeqs){
    int i;
    for(i = 0; i < nSeqs; i++){
        free(seqs[i].gid);
        free(seqs[i].str);
    } 

    return EXIT_SUCCESS;
}
