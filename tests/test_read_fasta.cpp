/**
 * @author jeff.daily@pnnl.gov
 *
 * Test the function read_fasta in libpgraph.
 */
#include "config.h"

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

#include <mpi.h>

#if HAVE_ARMCI
#include <mpi.h>
extern "C" {
#include <armci.h>
}
#endif

#if HAVE_ARMCI
#include "SequenceDatabaseArmci.hpp"
#else
#include "SequenceDatabaseReplicated.hpp"
#endif
#include "mpix.hpp"
#include "csequence.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exit;
using std::istringstream;
using std::vector;

using namespace pgraph;

static void sync_exit(int value, MPI_Comm comm)
{
        MPI_Comm_free(&comm);
        MPI_Finalize();
        exit(value);
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0;
    int nprocs = 0;
    long budget = 0;
    char budget_multiplier = 0;

    /* initialize MPI */
#if defined(THREADED)
    int provided = 0;
    MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided));
    assert(provided == MPI_THREAD_MULTIPLE);
#else
    MPI_CHECK(MPI_Init(&argc, &argv));
#endif
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));
    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &nprocs));

    /* sanity check that we got the correct number of arguments */
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

    {
#if HAVE_ARMCI
        SequenceDatabaseArmci sd(argv[1], budget, MPI_COMM_WORLD);
#else
        SequenceDatabaseReplicated sd(argv[1], budget, MPI_COMM_WORLD);
#endif
        sequences_t *sequences = pg_load_fasta(argv[1], '$');

        for (size_t i=0; i<sd.get_global_count(); ++i) {
            Sequence &seq1 = sd.get_sequence(i);
            sequence_t *seq2 = &sequences->seq[i];
            string str1(seq1);
            string str2(seq2->str);
            if (str1 != str2) {
                cout << "[" << rank << "] i=" << i << " mismatch" << endl;
                cout << str1 << endl;
                cout << str2 << endl;

            }
            assert(str1==str2);
        }

#if HAVE_ARMCI
        ARMCI_Barrier();
#else
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        pg_free_sequences(sequences);
    }

    /* clean up */
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
