/* C++ STL */
#include <string>
#include <vector>

/* 3rd party headers */
#include <mpi.h>

/* my headers */
#include "io.h"
#include "mpix.h"

using std::string;
using std::vector;


void read_fasta(const char *filename, size_t budget, vector<string> &sequences)
{
    MPI_Offset filesize;
    MPI_Offset localsize;
    MPI_Offset start;
    MPI_Offset end;
    MPI_Offset overlap=1024;
    MPI_File in;
    int rank;
    int size;
    char *chunk;
    int ierr;
    int i;

    sequences.clear();
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_CHECK(ierr);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_CHECK(ierr);
    ierr = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_CHECK(ierr);

    /* figure out who reads what */

    ierr = MPI_File_get_size(in, &filesize);
    MPI_CHECK(ierr);
    if (budget > filesize) {
        start = 0;
        end = filesize;
        localsize = filesize;
    }
    else {
        localsize = filesize/size;
        if (localsize > budget) {
            fprintf(stderr, "sequence memory budget not sufficient\n");
            MPI_Abort(MPI_COMM_WORLD, localsize-budget);
        }
        start = rank * localsize;
        /* we fudge the margins based on the memory budget specified */
        start = start - ((budget-localsize)/2);
        end   = end   + ((budget-localsize)/2);
    }

    /* except the last processor, of course */
    if (rank == size-1) end = filesize;
    /* except ranks near the front */
    if (start < 0) start = 0;

    localsize = end - start + 1;

    /* allocate memory */
    chunk = static_cast<char*>(malloc((localsize+1)*sizeof(char)));

    /* everyone reads in their part */
    MPI_File_read_at_all(in, start, chunk, localsize, MPI_CHAR,
            MPI_STATUS_IGNORE);
    chunk[localsize] = '\0';

    for (i=0; i<size; ++i) {
        if (rank == i) {
            printf("[%d] %s\n", rank, chunk);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    /*
     * everyone calculate what their start and end *really* are by going 
     * from the first newline after start to the first newline after the
     * overlap region starts (eg, after end - overlap + 1)
     */

#if 0
    int locstart=0, locend=localsize;
    if (rank != 0) {
        while(chunk[locstart] != '\n') locstart++;
        locstart++;
    }
    if (rank != size-1) {
        locend-=overlap;
        while(chunk[locend] != '\n') locend++;
    }
    localsize = locend-locstart+1;

    /* Now let's copy our actual data over into a new array, with no overlaps */
    char *data = (char *)malloc((localsize+1)*sizeof(char));
    memcpy(data, &(chunk[locstart]), localsize);
    data[localsize] = '\0';
    free(chunk);

    /* Now we'll count the number of lines */
    *nlines = 0;
    for (int i=0; i<localsize; i++)
        if (data[i] == '\n') (*nlines)++;

    /* Now the array lines will point into the data array at the start of each
     * line assuming nlines > 1 */
    *lines = (char **)malloc((*nlines)*sizeof(char *));
    (*lines)[0] = strtok(data,"\n");
    for (int i=1; i<(*nlines); i++)
        (*lines)[i] = strtok(NULL, "\n");
#endif

    return;
}
