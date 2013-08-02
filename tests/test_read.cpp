/**
 * @author jeff.daily@pnnl.gov
 *
 * Reading large data files collectively can be difficult. We don't want every
 * process opening and reading the file so we would like to have MPIIO handle
 * the details. But MPIIO can fail on systems with limited memory.
 */
#include "config.h"

#include <mpi.h>
#include <sys/stat.h>

#include <cassert>
#include <climits>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define ARG_LEN_MAX 1024

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

int rank = 0;
int nprocs = 0;
int check_count = 0;


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int provided;
    long file_size = -1;
    char *file_buffer = NULL;
    MPI_File fh;
    MPI_Status status;
    unsigned long seg_count=0;

    check_count = 0;

    /* initialize MPI */
#if defined(THREADED)
    MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided));
    assert(provided == MPI_THREAD_MULTIPLE);
#else
    MPI_CHECK(MPI_Init(&argc, &argv));
#endif
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));
    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &nprocs));

    /* sanity check that we got the correct number of arguments */
    if (argc <= 1 || argc >= 3) {
        if (0 == rank) {
            if (argc <= 1) {
                printf("missing input file\n");
            }
            else if (argc >= 3) {
                printf("too many arguments\n");
            }
            printf("usage: test_read sequence_file\n");
        }
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }

    /* we don't want all procs to stat the file! */
    if (0 == rank) {
        struct stat st;
        assert(0 == stat(argv[1], &st));
        file_size = st.st_size;
    }

    /* the file_size is broadcast to all */
    MPI_CHECK(MPI_Bcast(&file_size, 1, MPI_LONG, 0, comm));
    if (0 == rank) {
        printf("file_size=%ld\n", file_size);
    }

    if (0 == rank) {
        FILE *file = NULL;
        size_t read_count = 0;
        unsigned long my_seg_count = 0;
        char *file_buffer_dup = new char[file_size];

        file = fopen(argv[1], "r");
        if (NULL == file) {
            perror("fopen");
            printf("unable to open file on process 0\n");
            MPI_Abort(comm, 1);
        }
        printf("process 0 opened file\n");

        read_count = fread(file_buffer, file_size, 1, file);
        if (0 == read_count) {
            printf("unable to read file on process 0\n");
            MPI_Abort(comm, 1);
        }
        printf("process 0 read file\n");

        fclose(file);
        
        assert(file_buffer[0] == '>');
        ++my_seg_count;
        for (int i=1; i<file_size; ++i) {
            if (file_buffer[i] == '>' && file_buffer[i-1] == '\n') {
                ++my_seg_count;
            }
        }

        printf("[0] my_seg_count=%lu\n", my_seg_count);

        delete [] file_buffer_dup;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* allocate a buffer for the file, of the entire size */
    file_buffer = new char[file_size];
    char *current_file_buffer = file_buffer;
    MPI_Offset offset = 0;
    long file_size_remaining = file_size;

    /* all procs read the entire file */
    /* read the file in INT_MAX chunks due to MPI 'int' interface */
    MPI_CHECK(MPI_File_open(comm, argv[1],
                MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh));
    if (0 == rank) {
        cout << "INT_MAX=" << INT_MAX << endl;
    }
    //int chunk_size = INT_MAX;
    //int chunk_size = 8192;
    int chunk_size = 536870912; // 0.5 GB
    do {
        int count;
        int count_result;
        if (file_size_remaining > chunk_size) {
            count = chunk_size;
        }
        else {
            count = file_size_remaining;
        }
        MPI_CHECK(MPI_File_read_at_all(fh, offset, current_file_buffer,
                    count, MPI_CHAR, &status));
        MPI_CHECK(MPI_Get_count(&status, MPI_CHAR, &count_result));
        for (int p=0; p<nprocs; ++p) {
            if (p == rank) {
                cout << p << " read " << count_result << " bytes" << endl;
            }
            MPI_Barrier(comm);
        }
        offset += count_result;
        current_file_buffer += count_result;
        file_size_remaining -= count_result;
    }
    while (file_size_remaining > 0);
    MPI_CHECK(MPI_File_close(&fh));

    /* each process counts how many '>' characters are in the file_buffer */
    assert(file_buffer[0] == '>');
    ++seg_count;
    for (int i=1; i<file_size; ++i) {
        if (file_buffer[i] == '>' && file_buffer[i-1] == '\n') {
            ++seg_count;
#if 1
            if (seg_count % 10000 == 0) {
                MPI_Barrier(MPI_COMM_WORLD);
            }
#endif
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* print the seg_count on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] seg_count=%ld\n", rank, seg_count);
        }
        MPI_Barrier(comm);
    }



    /* clean up */
    delete [] file_buffer;
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
