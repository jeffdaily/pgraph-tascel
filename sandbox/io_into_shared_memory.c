/**
 * @author jeff.daily@pnnl.gov
 *
 * A first attempt at using MPI IO, reading on process 0 and storing the file
 * into a shared memory segment. It mostly worked, but I was advised by Sriram
 * to instead forget about master processes on each node, shared memory, and
 * the like.  Instead, I should use one MPI process per node and a threaded
 * model such that each thread has access to the same data structures.  Good
 * idea.
 *
 * DEPRECATED! But kept for reference, just in case...
 */
#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>


#define ARG_LEN_MAX 1024

int rank;
int nprocs;
int check_count;

#define MPI_CHECK(what) do {                              \
    int __err;                                            \
    __err = what;                                         \
    ++check_count;                                        \
    if (MPI_SUCCESS != __err) {                           \
        printf("[%d] FAILED FILE=%s LINE=%d:" #what "\n", \
                rank, __FILE__, __LINE__);                \
        MPI_Abort(comm, check_count);                     \
    }                                                     \
} while (0)


static void free_command_line(int argc, char **argv)
{
    int i;
    for (i=0; i<argc; ++i) {
        free(argv[i]);
        argv[i] = NULL;
    }
    free(argv);
    argv = NULL;
}


int main(int argc, char **argv)
{
    int all_argc = 0;
    char **all_argv = NULL;
    char *arg = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int i = 0;
    long file_size = -1;
    char *file_buffer = NULL;
    char *file_curbuf = NULL;
    MPI_File fh;
    MPI_Status status;
    int shmid = -1;
    long seg_count = 0;
    long *seg_idx;
    long *seg_len;
    int seg_label_found = 0;

    check_count = 0;

    MPI_CHECK(MPI_Init(&argc, &argv));
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));

    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank) {
        all_argc = argc;
    }
    MPI_CHECK(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm));
    all_argv = (char**)malloc(argc * sizeof(char*));
    for (i=0; i<all_argc; ++i) {
        int count = 0;
        if (0 == rank) {
            count = strlen(argv[i]) + 1;
        }
        MPI_CHECK(MPI_Bcast(&count, 1, MPI_INT, 0, comm));
        all_argv[i] = (char*)malloc(count * sizeof(char));
        if (0 == rank) {
            (void)strncpy(all_argv[i], argv[i], count);
        }
        MPI_CHECK(MPI_Bcast(all_argv[i], count, MPI_CHAR, 0, comm));
    }

#if DEBUG
    /* print the command line arguments */
    for (i=0; i<nprocs; ++i) {
        if (i == rank) {
            int j;
            for (j=0; j<all_argc; ++j) {
                printf("[%d] argv[%d]=%s\n", rank, j, all_argv[j]);
            }
        }
        MPI_Barrier(comm);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argc <= 1 || all_argc >= 3) {
        if (0 == rank) {
            printf("missing input file\n");
        }
        free_command_line(all_argc, all_argv);
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }

    /* process 0 open the file locally to determine its size */
    if (0 == rank) {
        FILE *file = fopen(all_argv[1], "r");
        if (NULL == file) {
            printf("unable to open file on process 0\n");
            MPI_Abort(comm, 1);
        }
        if (0 != fseek(file, 0, SEEK_END)) {
            printf("unable to seek to end of file on process 0\n");
            MPI_Abort(comm, 1);
        }
        file_size = ftell(file);
        if (-1 == file_size) {
            printf("unable to get size of file on process 0\n");
            MPI_Abort(comm, 1);
        }
        if (0 != fclose(file)) {
            printf("unable to get close file on process 0\n");
            MPI_Abort(comm, 1);
        }
    }
    /* process 0 (TODO: on each SMP node) creates a shared memory segment for
     * the input, then broadcasts within the local SMP node (TODO) */
    if (0 == rank) {
        errno = 0;
        shmid = shmget(IPC_PRIVATE, file_size, SHM_R|SHM_W);
        if (-1 == shmid) {
            perror("shmget");
            MPI_Abort(comm, -1);
        }
    }
    MPI_CHECK(MPI_Bcast(&shmid, 1, MPI_INT, 0, comm));
    MPI_Barrier(comm);
    errno = 0;
    if (0 == rank) {
        file_buffer = (char*)shmat(shmid, NULL, 0);
    }
    else {
        file_buffer = (char*)shmat(shmid, NULL, SHM_RDONLY);
    }
    if ((void*)-1 == file_buffer) {
        perror("shmat");
        MPI_Abort(comm, -1);
    }

    /* the file_size is broadcast to all, even though only master on each node
     * will need it
     * TODO: consider creating an MPI_Comm for all masters */
    MPI_CHECK(MPI_Bcast(&file_size, 1, MPI_LONG, 0, comm));
    MPI_CHECK(MPI_File_open(comm, all_argv[1],
                MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN,
                MPI_INFO_NULL, &fh));
    /* only process 0 reads the file, but all procs open it
     * TODO: See above about shared memory; also, the master process on each
     * node should be creating this memory and sharing it within the node */
    if (0 == rank) {
        printf("[%d] file_size=%ld\n", rank, file_size);
        MPI_CHECK(MPI_File_read_all(fh, file_buffer, file_size, MPI_CHAR, &status));
    }
    else {
        MPI_CHECK(MPI_File_read_all(fh, NULL, 0, MPI_CHAR, &status));
    }
    MPI_CHECK(MPI_File_close(&fh));

    /* each process counts how many '>' characters are in the file_buffer */
    for (i=0; i<file_size; ++i) {
        if (file_buffer[i] == '>') {
            ++seg_count;
        }
    }
#if 1
    /* print the seg_count on each process */
    for (i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] seg_count=%ld\n", rank, seg_count);
        }
        MPI_Barrier(comm);
    }
#endif
    seg_count = 0;
    seg_idx = (long*)malloc(seg_count * sizeof(long));
    seg_len = (long*)malloc(seg_count * sizeof(long));

#if 0
    /* each process indexes the file_buffer */
    i = 0;
    while (i<file_size) {
        if (file_buffer[i] == '>') {
            /* increment i until end of segment label line */
            while (file_buffer[i] != '\n') ++i;
            ++i; /* skip past the newline */
        }
        seg_idx[seg_count] = i;
        /* increment i until end of segment */
        while (file_buffer[i] != '\n') ++i;
        seg_idx[seg_count] = i - seg_idx[seg_count];
        ++seg_count;
        ++i;
    }
#endif

    /* clean up */
    errno = 0;
    if (-1 == shmdt(file_buffer)) {
        perror("shmdt");
    }
    if (0 == rank) {
        errno = 0;
        if (-1 == shmctl(shmid, IPC_RMID, 0)) {
            perror("shmctl(shmid, IPC_RMID, 0)");
        }
    }
    free(seg_idx);
    free(seg_len);
    free_command_line(all_argc, all_argv);
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
