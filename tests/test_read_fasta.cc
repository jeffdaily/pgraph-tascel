/**
 * @author jeff.daily@pnnl.gov
 *
 * Test the function read_fasta in libpgraph.
 */
#include <cstdlib>
#include <cassert>
#include <vector>

#include <mpi.h>

#include "SequenceDatabase.h"
#include "mpix.h"

using std::vector;


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int provided;
    int rank = 0;
    int nprocs = 0;
    vector<string> sequences;
    long budget;
    size_t budget_strlen;

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
    printf("argc=%d\n", argc);
    if (argc <= 2 || argc >= 4) {
        if (0 == rank) {
            if (argc <= 1) {
                printf("missing input file\n");
            }
            if (argc <= 2) {
                printf("missing memory budget\n");
            }
            else if (argc >= 4) {
                printf("too many arguments\n");
            }
            printf("usage: test_read_fasta sequence_file memory_budget\n");
        }
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }
    budget_strlen = strlen(argv[2]);
    printf("budget_strlen=%d\n", (int)budget_strlen);
    if (argv[2][budget_strlen-1] == 'm' || argv[2][budget_strlen-1] == 'M') {
        argv[2][budget_strlen-1]= '\0';
        budget = atol(argv[2]) * 1048576;
    }
    else if (argv[2][budget_strlen-1] == 'g' || argv[2][budget_strlen-1] == 'G') {
        argv[2][budget_strlen-1]= '\0';
        budget = atol(argv[2]) * 1073741824;
    }
    else {
        budget = atol(argv[2]);
    }
    printf("memory budget=%ld\n", budget);

    SequenceDatabase sd(argv[1], budget);

    /* clean up */
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
