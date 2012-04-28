/**
 * @author jeff.daily@pnnl.gov
 *
 * Reading large data files collectively can be difficult. We don't want every
 * process opening and reading the file so we would like to have MPIIO handle
 * the details. But MPIIO can fail on systems with limited memory.
 */
#include <mpi.h>
#include <sys/stat.h>

#include <cassert>
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
    vector<string> all_argv;
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

    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank) {
        MPI_CHECK(MPI_Bcast(&argc, 1, MPI_INT, 0, comm));
        for (int i=0; i<argc; ++i) {
            int length = strlen(argv[i])+1;
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(argv[i], length, MPI_CHAR, 0, comm));
            all_argv.push_back(argv[i]);
        }
    } else {
        int all_argc;
        MPI_CHECK(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm));
        for (int i=0; i<all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm));
            all_argv.push_back(buffer);
        }
    }

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 3) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() >= 3) {
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
        assert(0 == stat(all_argv[0].c_str(), &st));
        file_size = st.st_size;
    }

    /* the file_size is broadcast to all */
    MPI_CHECK(MPI_Bcast(&file_size, 1, MPI_LONG, 0, comm));
    if (0 == rank) {
        printf("file_size=%ld\n", file_size);
    }

    /* allocate a buffer for the file, of the entire size */
    file_buffer = new char[file_size];
    MPI_Datatype newtype = MPI_CHAR;
    int newsize = 0;
    if (file_size > sizeof(int)) {
        if (0 == file_size % sizeof(long long)) {
            newtype = MPI_LONG_LONG_INT;
            newsize = file_size / sizeof(long long);
        }
        else if (0 == file_size % sizeof(long)) {
            newtype = MPI_LONG;
            newsize = file_size / sizeof(long);
        }
        else if (0 == file_size % sizeof(int)) {
            newtype = MPI_INT;
            newsize = file_size / sizeof(int);
        }
        else if (0 == file_size % sizeof(short)) {
            newtype = MPI_SHORT;
            newsize = file_size / sizeof(short);
        }
        else {
            if (0 == rank) {
                cout << "odd file size > 2GB " << file_size << endl;
            }
            MPI_Finalize();
            return 1;
        }
    }
    if (newsize > sizeof(int)) {
        if (0 == rank) {
            cout << "file size still > 2GB " << newsize << endl;
        }
        MPI_Finalize();
        return 1;
    }

    /* all procs read the entire file */
    MPI_CHECK(MPI_File_open(comm, const_cast<char*>(all_argv[1].c_str()),
                MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh));
    MPI_CHECK(MPI_File_read_all(fh, file_buffer, newsize, newtype, &status));
    MPI_CHECK(MPI_File_close(&fh));

    /* each process counts how many '>' characters are in the file_buffer */
    for (int i=0; i<file_size; ++i) {
        if (file_buffer[i] == '>') {
            ++seg_count;
        }
    }

    /* print the seg_count on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] seg_count=%ld\n", rank, seg_count);
        }
        MPI_Barrier(comm);
    }

    /* clean up */
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
