#ifndef MPIX_H_
#define MPIX_H_

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;

extern int check_count;

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

void mpix_bcast_argv(MPI_Comm comm, int argc, char **argv, vector<string> &all_argv);
void mpix_read_file(MPI_Comm comm, const string &file_name, char* &file_buffer, long &file_size);

void mpix_print_sync(MPI_Comm comm, const string &name, const vector<string> &what);
void mpix_print_sync(MPI_Comm comm, const string &name, const string &what);

template <class T>
void mpix_print_sync(MPI_Comm comm, const string &name, const T &what)
{
    int rank;
    int size;
    T what_copy = what;
    T *all_what = NULL;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &size));

    all_what = new T[size];
    MPI_CHECK(MPI_Gather(&what_copy, sizeof(T), MPI_CHAR, all_what, sizeof(T), MPI_CHAR, 0, comm));
    if (0 == rank) {
        for (int i=0; i<size; ++i) {
            cout << "[" << i << "] " << name << "=" << all_what[i] << endl;
        }
    }
    delete [] all_what;
    MPI_Barrier(comm);
}

#endif /* MPIX_H_ */
