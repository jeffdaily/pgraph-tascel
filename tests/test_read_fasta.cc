/**
 * @author jeff.daily@pnnl.gov
 *
 * Test the function read_fasta in libpgraph.
 */
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

#include <mpi.h>

#include "SequenceDatabase.h"
#include "mpix.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exit;
using std::istringstream;
using std::vector;


static void sync_exit(int value, MPI_Comm comm)
{
        MPI_Comm_free(&comm);
        MPI_Finalize();
        exit(value);
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int provided = 0;
    int rank = 0;
    int nprocs = 0;
    vector<string> sequences;
    long budget = 0;
    char budget_multiplier = 0;

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
    cout << "argc=" << argc << endl;
    if (argc <= 2 || argc >= 4) {
        if (0 == rank) {
            if (argc <= 1) {
                cerr << "missing input file" << endl;
            }
            if (argc <= 2) {
                cerr << "missing memory budget" << endl;
            }
            else if (argc >= 4) {
                cerr << "too many arguments" << endl;
            }
            cerr << "usage: " << argv[0]
                << " sequence_file memory_budget" << endl;
        }
        sync_exit(1, comm);
    }
    
    {
        istringstream iss(argv[2]);
        iss >> budget >> budget_multiplier;
        if (budget <= 0) {
            cerr << "memory budget must be positive real number" << endl;
            sync_exit(1, comm);
        }
    }

    if (budget_multiplier == 'k' || budget_multiplier == 'K') {
        budget *= 1024; /* kilobyte */
    }
    else if (budget_multiplier == 'm' || budget_multiplier == 'M') {
        budget *= 1048576; /* megabyte */
    }
    else if (budget_multiplier == 'g' || budget_multiplier == 'G') {
        budget *= 1073741824; /* gigabyte */
    }
    cout << "memory budget=" << budget << " bytes" << endl;

    SequenceDatabase sd(argv[1], budget);
    mpix_print_sync(comm, "local_count", sd.get_local_count());
    mpix_print_sync(comm, "global_count", sd.get_global_count());
    Sequence &s1 = sd.get_sequence(0);
    mpix_print_sync(comm, "seq 0", string(s1));
    Sequence &s2 = sd.get_sequence(1);
    mpix_print_sync(comm, "seq 0", string(s2));
    int score,ndig,alen;
    s1.align(s2, score, ndig, alen);
    mpix_print_sync(comm, "score", score);
    mpix_print_sync(comm, "ndig", ndig);
    mpix_print_sync(comm, "alen", alen);

    /* clean up */
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
