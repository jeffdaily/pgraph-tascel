#include "loadseq.h"
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
 *------------------------------------------------*/
void loadAllSeqs(char *fileName, int nSeqs, SEQ *seqs, long *nChars){
    char *buffer;
    struct stat buf;
    luc_error_t err;
    int snap_err = 0;

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
    int *flags = malloc(nSeqs * sizeof(int));

#pragma mta assert no dependence
    for (int i = 0; i < buffer_size; i++)
        if (buffer[i] == FASTA_FLAG) flags[ int_fetch_add(&next_flag, 1) ] = i;

#pragma mta assert no dependence
    for (int i = 0; i < nSeqs; i++) {
        
        int ndx = flags[i] + 1;
        seqs[i].gid = buffer + ndx;

        while (! isspace(buffer[ndx])) ndx++;
        buffer[ndx] = '\0';

        ndx++;
        int strLen  = 0;
        seqs[i].str = buffer + ndx;
        while (! isspace(buffer[ndx])) {ndx++; strLen++;}

        buffer[ndx]     = DOLLAR;
        buffer[ndx + 1] = '\0';
        seqs[i].strLen  = strLen;
        NN             += strLen;
    }

    *nChars = NN;
}

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
