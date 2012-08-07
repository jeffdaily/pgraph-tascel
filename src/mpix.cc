#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "mpix.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int check_count=0;


/* MPI standard does not guarantee all procs receive argc and argv */
void mpix_bcast_argv(
        MPI_Comm comm, int argc, char **argv, vector<string> &all_argv)
{
    int rank;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));

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
}


void mpix_print_sync(
        MPI_Comm comm, const string &name, const vector<string> &what)
{
    int rank;
    int size;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &size));

    for (int i=0; i<size; ++i) {
        if (i == rank) {
            int j;
            for (j=0; j<what.size(); ++j) {
                cout << "[" << rank << "] "
                    << name << "[" << j << "]="
                    << what.at(j) << endl;
            }
        }
        MPI_Barrier(comm);
    }
}


void mpix_print_sync(MPI_Comm comm, const string &name, const string &what)
{
    int rank;
    int size;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &size));

    for (int i=0; i<size; ++i) {
        if (i == rank) {
            cout << "[" << rank << "] " << name << "[" << i << "]=" << what << endl;
        }
        MPI_Barrier(comm);
    }
}


/**
 * Collectively read a file.
 *
 * All processes within the MPI_Comm instance will receive the entire contents
 * of the file.
 *
 * @param[in] comm instance
 * @param[in] file_name to open
 * @param[out] file_buffer to store file contents
 * @param[out] file_size of the given file
 */
void mpix_read_file(
        MPI_Comm comm, const string &file_name,
        char* &file_buffer, long &file_size)
{
    int rank;
    int size;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &size));

    /* process 0 stats the file to determine its size */
    if (0 == rank) {
        struct stat statbuf;
        int retval;

        retval = stat(file_name.c_str(), &statbuf);
        if (-1 == retval) {
            perror("stat");
            printf("unable to stat sequence file\n");
            MPI_Abort(comm, 1);
        }
        file_size = statbuf.st_size;
    }

    /* the file_size is broadcast to all */
    MPI_CHECK(MPI_Bcast(&file_size, 1, MPI_LONG, 0, comm));
    if (0 == rank) {
        printf("file_size=%ld\n", file_size);
    }

    /* allocate a buffer for the file, of the entire size */
    file_buffer = new char[file_size];

#if MPIX_READ_FILE_ALL
    /* all procs read the entire file */
    MPI_CHECK(MPI_File_open(comm, const_cast<char*>(file_name.c_str()),
                MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN,
                MPI_INFO_NULL, &fh));
    MPI_CHECK(MPI_File_read_all(fh, file_buffer, file_size, MPI_CHAR, &status));
    MPI_CHECK(MPI_File_close(&fh));
#else
    if (0 == rank) {
        FILE *file = NULL;
        size_t read_count = 0;

        file = fopen(file_name.c_str(), "r");
        if (NULL == file) {
            perror("fopen");
            printf("unable to open file on process 0\n");
            MPI_Abort(comm, 1);
        }

        read_count = fread(file_buffer, file_size, 1, file);
        if (0 == read_count) {
            printf("unable to read file on process 0\n");
            MPI_Abort(comm, 1);
        }

    }

#   define BCAST_IN_CHUNKS 1
#   if BCAST_IN_CHUNKS
    /* bcast file contents in 1GB chunks */
    long chunk_size = 1073741824;
    long offset = 0;
    while (offset < file_size) {
        long message_size = chunk_size;
        if (offset+chunk_size > file_size) {
            message_size = file_size % chunk_size;
        }
        MPI_CHECK(MPI_Bcast(&file_buffer[offset],
                    message_size, MPI_CHAR, 0, comm));
        offset += chunk_size;
    }
#   else
    /* bcast file contents */
    MPI_CHECK(MPI_Bcast(file_buffer, file_size, MPI_CHAR, 0, comm));
#   endif
#endif
}
