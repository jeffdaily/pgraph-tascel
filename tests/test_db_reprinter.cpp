/**
 * @author jeff.daily@pnnl.gov
 *
 * Reads in fasta file and reprints. This test was designed for diffing the
 * resulting file agains the original in cases where a distributed DB is used.
 */
#include "config.h"

#include <mpi.h>
#include <tascel.h>
#if THREADED
#include "tascelx.hpp"
#endif

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "Bootstrap.hpp"
#include "mpix.hpp"
#include "mpix_types.hpp"
#include "Parameters.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"
#include "SequenceDatabaseTascel.hpp"
#include "SequenceDatabaseWithStats.hpp"

using namespace std;
using namespace tascel;
using namespace pgraph;

int main(int argc, char **argv)
{
    double time_main = 0.0;
    int rank = 0;
    int nprocs = 0;
    vector<string> all_argv;
    SequenceDatabase *db = NULL;
    Parameters *parameters = NULL;
    char delimiter = '\0';

    pgraph::initialize(argc, argv);
    rank = mpix::comm_rank(pgraph::comm);
    nprocs = mpix::comm_size(pgraph::comm);
    time_main = MPI_Wtime();

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, pgraph::comm);
#if THREADED
    allocate_threads();
    initialize_server_thread();
#endif
    
    parameters = new Parameters;

    all_argv = mpix::bcast(argc, argv, pgraph::comm);

    /* print the command line arguments */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "--- Command Line ";
        cout << header.str() << endl;
        for (size_t j=0; j<all_argv.size(); ++j) {
            cout << "argv[" << j << "]='" << all_argv[j] << "'" << endl;
        }
    }

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                cout << "missing input file" << endl;
            }
            else if (all_argv.size() >= 4) {
                cout << "too many arguments" << endl;
            }
            cout << "usage: align sequence_file <config_file>" << endl;
        }
        TascelConfig::finalize();
        pgraph::finalize();
        return 1;
    }
    else if (all_argv.size() >= 3) {
        parameters->parse(all_argv[2].c_str(), pgraph::comm);
    }
    else if (all_argv.size() >= 2) {
        /* do nothing */
    }

    /* print parameters */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "-- Parameters ";
        cout << header.str() << endl;
        cout << *parameters << endl;
    }

    if (parameters->use_tree
            || parameters->use_tree_dynamic
            || parameters->use_tree_hybrid) {
        delimiter = parameters->alphabet_dollar; /* dollar */
    }

    if (0 == rank) {
        cout << "opening file " << all_argv[1] << endl;
    }
    db = new SequenceDatabaseTascel(all_argv[1],
            parameters->memory_sequences, pgraph::comm, delimiter);
    if (0 == rank) {
        cout << "opened file " << all_argv[1] << endl;
    }

#if THREADED
    if (0 == rank) {
        cout << "allocating threads" << endl;
    }
    allocate_threads();
    if (0 == rank) {
        cout << "initializing server thread" << endl;
    }
    initialize_server_thread();
#endif

    if (0 == rank) {
        ofstream out("reprint.out");
        for (size_t i=0; i<db->size(); ++i) {
            string id;
            string s;
            Sequence *sequence = NULL;
            
            cout << "getting sequence " << i << endl;
            sequence = db->get_sequence(i);
            sequence->get_id(id);
            sequence->get_sequence(s);
            cout << "writing sequence " << i << endl;
            out << id << endl;
            if (sequence->uses_delimiter()) {
                out << s.substr(0,s.length()-1) << endl;
            }
            else {
                out << s << endl;
            }
            delete sequence;
        }
        cout << "closing file" << endl;
        out.close();
    }

#if THREADED
    if (0 == rank) {
        cout << "finalizing server thread" << endl;
    }
    finalize_server_thread();
    if (0 == rank) {
        cout << "deallocating threads" << endl;
    }
    deallocate_threads();
#endif

#if 0
    if (0 == rank) {
        cout << "amBarrier" << endl;
    }
    amBarrier();
    if (0 == rank) {
        cout << "MPI_Barrier" << endl;
    }
    MPI_Barrier(pgraph::comm);
#endif

    time_main = MPI_Wtime() - time_main;
    if (0 == rank) {
        cout << "time_main " << time_main << " seconds" << endl;
    }
    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}
