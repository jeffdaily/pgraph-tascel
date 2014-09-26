#include "config.h"

#include <stdint.h>

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "align/align.h"
#include "blosum/blosum62.h"
#include "Bootstrap.hpp"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"

using namespace ::std;
using namespace ::pgraph;

#define ENABLE_WOZNIAK 1
#define ENABLE_STRIPED 1

static float pct(float orig_, float new_)
{
    return 100.0*(orig_-new_)/orig_;
}


int main(int argc, char **argv)
{
    int rank = 0;
    int nprocs = 0;
    SequenceDatabase *db = NULL;
    size_t db_size = 0U;
    size_t db_longest = 0U; /* +16 for woz padding */;
    Parameters *parameters = NULL;
    int ref_score = 0;
    int ref_matches = 0;
    int ref_length = 0;
    int score = 0;
    int matches = 0;
    int length = 0;
    size_t limit = 10000U;
    size_t i = 0U;
    int * tbl_pr = NULL;
    int * del_pr = NULL;
    int * mch_pr = NULL;
    int * len_pr = NULL;
    vector<string> all_argv;
    char delimiter = '\0';
    double time = 0.0f;
    double time_nw_ref = 0.0f;
    double time_nw_woz = 0.0f;
    double time_nw_striped = 0.0f;
    double time_sg_ref = 0.0f;
    double time_sg_woz = 0.0f;
    double time_sg_striped = 0.0f;
    double time_sw_ref = 0.0f;
    double time_sw_woz = 0.0f;
    double time_sw_striped = 0.0f;
    double time_nw_stats_ref = 0.0f;
    double time_nw_stats_woz = 0.0f;
    double time_nw_stats_striped = 0.0f;
    double time_sg_stats_ref = 0.0f;
    double time_sg_stats_woz = 0.0f;
    double time_sg_stats_striped = 0.0f;
    double time_sw_stats_ref = 0.0f;
    double time_sw_stats_woz = 0.0f;
    double time_sw_stats_striped = 0.0f;

    /* MPI standard does not guarantee all procs receive argc and argv */
    pgraph::initialize(argc, argv);
    rank = mpix::comm_rank(pgraph::comm);
    nprocs = mpix::comm_size(pgraph::comm);
    all_argv = mpix::bcast(argc, argv, pgraph::comm);

    assert(0 == rank);
    assert(1 == nprocs);

    parameters = new Parameters;

    /* print the command line arguments */
    {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "--- Command Line ";
        cout << header.str() << endl;
        for (size_t j=0; j<all_argv.size(); ++j) {
            cout << "argv[" << j << "]='" << all_argv[j] << "'" << endl;
        }
        cout << endl;
    }

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (all_argv.size() <= 1) {
            cout << "missing input file" << endl;
        }
        else if (all_argv.size() >= 4) {
            cout << "too many arguments" << endl;
        }
        cout << "usage: align sequence_file <config_file>" << endl;
        return 1;
    }
    else if (all_argv.size() >= 3) {
        parameters->parse(all_argv[2].c_str(), pgraph::comm);
    }
    else if (all_argv.size() >= 2) {
        /* do nothing */
    }

    /* print parameters */
    {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "-- Parameters ";
        cout << header.str() << endl;
        cout << *parameters << endl;
        cout << endl;
    }

    time = MPI_Wtime();
    db = new SequenceDatabaseReplicated(all_argv[1],
            parameters->memory_sequences, pgraph::comm, delimiter);
    time = MPI_Wtime() - time;
    cout << "time sequence db open " << time << endl << endl;
    db_longest = db->longest() + 16; /* +16 for woz padding */
    db_size = db->size();

    /* allocate buffers for alignments */
    tbl_pr = new int[db_longest];
    del_pr = new int[db_longest];
    mch_pr = new int[db_longest];
    len_pr = new int[db_longest];

    for (size_t i=0; i<db_size; ++i) {
        for (size_t j=i+1; j<db_size; ++j) {
            Sequence *seqi = db->get_sequence(i);
            Sequence *seqj = db->get_sequence(j);
            const char * seqa = NULL;
            const char * seqb = NULL;
            size_t lena = 0U;
            size_t lenb = 0U;
            seqi->get_sequence(seqa, lena);
            seqj->get_sequence(seqb, lenb);

            time = MPI_Wtime();
            ref_score = nw(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_nw_ref += time;

#if ENABLE_WOZNIAK
            time = MPI_Wtime();
            score = nw_woz(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_nw_woz += time;

            assert(ref_score == score);
            if (ref_score != score) {
                cout << "nw woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED
            time = MPI_Wtime();
            score = nw_striped(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_nw_striped += time;

            assert(ref_score == score);
            if (ref_score != score) {
                cout << "nw stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = nw_stats(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &ref_matches, &ref_length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_nw_stats_ref += time;

#if ENABLE_WOZNIAK
            time = MPI_Wtime();
            score = nw_stats_woz(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_nw_stats_woz += time;

            assert(ref_score == score);
            if (ref_score != score
                    || ref_matches != matches
                    || ref_length != length) {
                cout << "nw stats woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED
            time = MPI_Wtime();
            score = nw_stats_striped(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_nw_stats_striped += time;

            assert(ref_score == score);
            if (ref_score != score
                    || ref_matches != matches
                    || ref_length != length) {
                cout << "nw stats stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = sg(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_sg_ref += time;

#if ENABLE_WOZNIAK
            time = MPI_Wtime();
            score = sg_woz(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_sg_woz += time;

            assert(ref_score == score);
            if (ref_score != score) {
                cout << "sg woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED
            time = MPI_Wtime();
            score = sg_striped(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sg_striped += time;

            assert(ref_score == score);
            if (ref_score != score) {
                cout << "sg stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = sg_stats(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &ref_matches, &ref_length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sg_stats_ref += time;

#if ENABLE_WOZNIAK
            time = MPI_Wtime();
            score = sg_stats_woz(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sg_stats_woz += time;

            assert(ref_score == score);
            if (ref_score != score
                    || ref_matches != matches
                    || ref_length != length) {
                cout << "sg stats woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED
            time = MPI_Wtime();
            score = sg_stats_striped(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sg_stats_striped += time;

            assert(ref_score == score);
            if (ref_score != score
                    || ref_matches != matches
                    || ref_length != length) {
                cout << "sg stats stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = sw(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_sw_ref += time;

#if ENABLE_WOZNIAK
            time = MPI_Wtime();
            score = sw_woz(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_sw_woz += time;

            assert(ref_score == score);
            if (ref_score != score) {
                cout << "sw woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED
            time = MPI_Wtime();
            score = sw_striped(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sw_striped += time;

            assert(ref_score == score);
            if (ref_score != score) {
                cout << "sw stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = sw_stats(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &ref_matches, &ref_length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sw_stats_ref += time;

#if ENABLE_WOZNIAK
            time = MPI_Wtime();
            score = sw_stats_woz(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sw_stats_woz += time;

            assert(ref_score == score);
            if (ref_score != score
                    || ref_matches != matches
                    || ref_length != length) {
                cout << "sw stats woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED
            time = MPI_Wtime();
            score = sw_stats_striped(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sw_stats_striped += time;

            assert(ref_score == score);
            if (ref_score != score
                    || ref_matches != matches
                    || ref_length != length) {
                cout << "sw stats stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

            delete seqi;
            delete seqj;
        }
    }

    cout << "nw reference\t" << time_nw_ref << endl;
#if ENABLE_WOZNIAK
    cout << "nw wozniak  \t" << time_nw_woz << "\t" << pct(time_nw_ref, time_nw_woz) << endl;
#endif
#if ENABLE_STRIPED
    cout << "nw striped  \t" << time_nw_striped << endl;
#endif

    cout << "sg reference\t" << time_sg_ref << endl;
#if ENABLE_WOZNIAK
    cout << "sg wozniak  \t" << time_sg_woz << "\t" << pct(time_sg_ref, time_sg_woz) << endl;
#endif
#if ENABLE_STRIPED
    cout << "sg striped  \t" << time_sg_striped << endl;
#endif

    cout << "sw reference\t" << time_sw_ref << endl;
#if ENABLE_WOZNIAK
    cout << "sw wozniak  \t" << time_sw_woz << "\t" << pct(time_sw_ref, time_sw_woz) << endl;
#endif
#if ENABLE_STRIPED
    cout << "sw striped  \t" << time_sw_striped << endl;
#endif

    cout << "nw stats reference\t" << time_nw_stats_ref << endl;
#if ENABLE_WOZNIAK
    cout << "nw stats wozniak  \t" << time_nw_stats_woz << "\t" << pct(time_nw_stats_ref, time_nw_stats_woz) << endl;
#endif
#if ENABLE_STRIPED
    cout << "nw stats striped  \t" << time_nw_stats_striped << endl;
#endif

    cout << "sg stats reference\t" << time_sg_stats_ref << endl;
#if ENABLE_WOZNIAK
    cout << "sg stats wozniak  \t" << time_sg_stats_woz << "\t" << pct(time_sg_stats_ref, time_sg_stats_woz) << endl;
#endif
#if ENABLE_STRIPED
    cout << "sg stats striped  \t" << time_sg_stats_striped << endl;
#endif

    cout << "sw stats reference\t" << time_sw_stats_ref << endl;
#if ENABLE_WOZNIAK
    cout << "sw stats wozniak  \t" << time_sw_stats_woz << "\t" << pct(time_sw_stats_ref, time_sw_stats_woz) << endl;
#endif
#if ENABLE_STRIPED
    cout << "sw stats striped  \t" << time_sw_stats_striped << endl;
#endif

    pgraph::finalize();

    delete [] tbl_pr;
    delete [] del_pr;
    delete [] mch_pr;
    delete [] len_pr;
    delete parameters;
    delete db;

    return 0;
}
