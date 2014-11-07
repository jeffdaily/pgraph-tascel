#include "config.h"

#include <mpi.h>

#include <stdint.h>

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "align/align.h"
#include "blosum/blosum62.h"
#include "Bootstrap.hpp"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"

#if HAVE_SSE2
#include "align/align_wozniak_128_16.h"
#include "align/align_striped_128_16.h"
#include "align/align_scan_128_16.h"
#endif

#if HAVE_SSE41
#include "align/align_scan_128_8.h"
#endif

using namespace ::std;
using namespace ::pgraph;

#define ENABLE_WOZNIAK 1
#define ENABLE_STRIPED 0
#define ENABLE_SCAN 1
#define STATS_MISMATCH 0
#define SATURATION 0

static float pct(float orig_, float new_)
{
    //return 100.0*(orig_-new_)/orig_;
    return orig_/new_;
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
    double time_nw_scan = 0.0f;
    double time_nw_scan8 = 0.0f;
    double time_sg_ref = 0.0f;
    double time_sg_woz = 0.0f;
    double time_sg_striped = 0.0f;
    double time_sg_scan = 0.0f;
    double time_sg_scan8 = 0.0f;
    double time_sw_ref = 0.0f;
    double time_sw_woz = 0.0f;
    double time_sw_striped = 0.0f;
    double time_sw_scan = 0.0f;
    double time_sw_scan8 = 0.0f;
    double time_nw_stats_ref = 0.0f;
    double time_nw_stats_woz = 0.0f;
    double time_nw_stats_striped = 0.0f;
    double time_nw_stats_scan = 0.0f;
    double time_nw_stats_scan8 = 0.0f;
    double time_sg_stats_ref = 0.0f;
    double time_sg_stats_woz = 0.0f;
    double time_sg_stats_striped = 0.0f;
    double time_sg_stats_scan = 0.0f;
    double time_sg_stats_scan8 = 0.0f;
    double time_sw_stats_ref = 0.0f;
    double time_sw_stats_woz = 0.0f;
    double time_sw_stats_striped = 0.0f;
    double time_sw_stats_scan = 0.0f;
    double time_sw_stats_scan8 = 0.0f;

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

#if ENABLE_WOZNIAK && HAVE_SSE2
            time = MPI_Wtime();
            score = nw_wozniak_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_nw_woz += time;

            if (ref_score != score) {
                cout << "nw woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED && HAVE_SSE2
            time = MPI_Wtime();
            score = nw_striped_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_nw_striped += time;

            if (ref_score != score) {
                cout << "nw stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE2
            time = MPI_Wtime();
            score = nw_scan_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_nw_scan += time;

            if (ref_score != score) {
                cout << "nw scan mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE41
            time = MPI_Wtime();
            score = nw_scan_128_8(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_nw_scan8 += time;

            if (SATURATION && ref_score != score) {
                cout << "nw scan8 mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = nw_stats(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &ref_matches, &ref_length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_nw_stats_ref += time;

#if ENABLE_WOZNIAK && HAVE_SSE2
            time = MPI_Wtime();
            score = nw_stats_wozniak_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_nw_stats_woz += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "nw stats woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED && HAVE_SSE2
            time = MPI_Wtime();
            score = nw_stats_striped_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_nw_stats_striped += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "nw stats stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE2
            time = MPI_Wtime();
            score = nw_stats_scan_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_nw_stats_scan += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "nw stats scan mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE41
            time = MPI_Wtime();
            score = nw_stats_scan_128_8(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_nw_stats_scan8 += time;

            if (SATURATION && ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "nw stats scan8 mismatch " << i << " " << j << ":\t";
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

#if ENABLE_WOZNIAK && HAVE_SSE2
            time = MPI_Wtime();
            score = sg_wozniak_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_sg_woz += time;

            if (ref_score != score) {
                cout << "sg woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED && HAVE_SSE2
            time = MPI_Wtime();
            score = sg_striped_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sg_striped += time;

            if (ref_score != score) {
                cout << "sg stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE2
            time = MPI_Wtime();
            score = sg_scan_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sg_scan += time;

            if (ref_score != score) {
                cout << "sg scan mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE41
            time = MPI_Wtime();
            score = sg_scan_128_8(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sg_scan8 += time;

            if (SATURATION && ref_score != score) {
                cout << "sg scan8 mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = sg_stats(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &ref_matches, &ref_length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sg_stats_ref += time;

#if ENABLE_WOZNIAK && HAVE_SSE2
            time = MPI_Wtime();
            score = sg_stats_wozniak_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sg_stats_woz += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sg stats woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED && HAVE_SSE2
            time = MPI_Wtime();
            score = sg_stats_striped_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sg_stats_striped += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sg stats stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE2
            time = MPI_Wtime();
            score = sg_stats_scan_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sg_stats_scan += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sg stats scan mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE41
            time = MPI_Wtime();
            score = sg_stats_scan_128_8(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sg_stats_scan8 += time;

            if (SATURATION && ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sg stats scan8 mismatch " << i << " " << j << ":\t";
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

#if ENABLE_WOZNIAK && HAVE_SSE2
            time = MPI_Wtime();
            score = sw_wozniak_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62, tbl_pr, del_pr);
            time = MPI_Wtime() - time;
            time_sw_woz += time;

            if (ref_score != score) {
                cout << "sw woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED && HAVE_SSE2
            time = MPI_Wtime();
            score = sw_striped_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sw_striped += time;

            if (ref_score != score) {
                cout << "sw stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE2
            time = MPI_Wtime();
            score = sw_scan_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sw_scan += time;

            if (ref_score != score) {
                cout << "sw scan mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE41
            time = MPI_Wtime();
            score = sw_scan_128_8(seqa, lena, seqb, lenb, 10, 1, blosum62__);
            time = MPI_Wtime() - time;
            time_sw_scan8 += time;

            if (SATURATION && ref_score != score) {
                cout << "sw scan8 mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << endl;
            }
#endif

            time = MPI_Wtime();
            ref_score = sw_stats(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &ref_matches, &ref_length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sw_stats_ref += time;

#if ENABLE_WOZNIAK && HAVE_SSE2
            time = MPI_Wtime();
            score = sw_stats_wozniak_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62,
                    &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
            time = MPI_Wtime() - time;
            time_sw_stats_woz += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sw stats woz mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_STRIPED && HAVE_SSE2
            time = MPI_Wtime();
            score = sw_stats_striped_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sw_stats_striped += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sw stats stp mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE2
            time = MPI_Wtime();
            score = sw_stats_scan_128_16(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sw_stats_scan += time;

            if (ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sw stats scan mismatch " << i << " " << j << ":\t";
                cout << " " << ref_score << "|" << score;
                cout << " " << ref_matches << "|" << matches;
                cout << " " << ref_length << "|" << length;
                cout << endl;
            }
#endif

#if ENABLE_SCAN && HAVE_SSE41
            time = MPI_Wtime();
            score = sw_stats_scan_128_8(seqa, lena, seqb, lenb, 10, 1, blosum62__,
                    &matches, &length);
            time = MPI_Wtime() - time;
            time_sw_stats_scan8 += time;

            if (SATURATION && ref_score != score
#if STATS_MISMATCH
                    || ref_matches != matches
                    || ref_length != length
#endif
                    ) {
                cout << "sw stats scan8 mismatch " << i << " " << j << ":\t";
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
#if ENABLE_WOZNIAK && HAVE_SSE2
    cout << "nw wozniak  \t" << time_nw_woz << "\t" << pct(time_nw_ref, time_nw_woz) << endl;
#endif
#if ENABLE_STRIPED && HAVE_SSE2
    cout << "nw striped  \t" << time_nw_striped << "\t" << pct(time_nw_ref, time_nw_striped) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE2
    cout << "nw scan     \t" << time_nw_scan << "\t" << pct(time_nw_ref, time_nw_scan) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE41
    cout << "nw scan8    \t" << time_nw_scan8 << "\t" << pct(time_nw_ref, time_nw_scan8) << endl;
#endif

    cout << "sg reference\t" << time_sg_ref << endl;
#if ENABLE_WOZNIAK && HAVE_SSE2
    cout << "sg wozniak  \t" << time_sg_woz << "\t" << pct(time_sg_ref, time_sg_woz) << endl;
#endif
#if ENABLE_STRIPED && HAVE_SSE2
    cout << "sg striped  \t" << time_sg_striped << "\t" << pct(time_sg_ref, time_sg_striped) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE2
    cout << "sg scan     \t" << time_sg_scan << "\t" << pct(time_sg_ref, time_sg_scan) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE41
    cout << "sg scan8    \t" << time_sg_scan8 << "\t" << pct(time_sg_ref, time_sg_scan8) << endl;
#endif

    cout << "sw reference\t" << time_sw_ref << endl;
#if ENABLE_WOZNIAK && HAVE_SSE2
    cout << "sw wozniak  \t" << time_sw_woz << "\t" << pct(time_sw_ref, time_sw_woz) << endl;
#endif
#if ENABLE_STRIPED && HAVE_SSE2
    cout << "sw striped  \t" << time_sw_striped << "\t" << pct(time_sw_ref, time_sw_striped) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE2
    cout << "sw scan     \t" << time_sw_scan << "\t" << pct(time_sw_ref, time_sw_scan) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE41
    cout << "sw scan8    \t" << time_sw_scan8 << "\t" << pct(time_sw_ref, time_sw_scan8) << endl;
#endif

    cout << "nw stats reference\t" << time_nw_stats_ref << endl;
#if ENABLE_WOZNIAK && HAVE_SSE2
    cout << "nw stats wozniak  \t" << time_nw_stats_woz << "\t" << pct(time_nw_stats_ref, time_nw_stats_woz) << endl;
#endif
#if ENABLE_STRIPED && HAVE_SSE2
    cout << "nw stats striped  \t" << time_nw_stats_striped << "\t" << pct(time_nw_stats_ref, time_nw_stats_striped) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE2
    cout << "nw stats scan     \t" << time_nw_stats_scan << "\t" << pct(time_nw_stats_ref, time_nw_stats_scan) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE41
    cout << "nw stats scan8    \t" << time_nw_stats_scan8 << "\t" << pct(time_nw_stats_ref, time_nw_stats_scan8) << endl;
#endif

    cout << "sg stats reference\t" << time_sg_stats_ref << endl;
#if ENABLE_WOZNIAK && HAVE_SSE2
    cout << "sg stats wozniak  \t" << time_sg_stats_woz << "\t" << pct(time_sg_stats_ref, time_sg_stats_woz) << endl;
#endif
#if ENABLE_STRIPED && HAVE_SSE2
    cout << "sg stats striped  \t" << time_sg_stats_striped << "\t" << pct(time_sg_stats_ref, time_sg_stats_striped) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE2
    cout << "sg stats scan     \t" << time_sg_stats_scan << "\t" << pct(time_sg_stats_ref, time_sg_stats_scan) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE41
    cout << "sg stats scan8    \t" << time_sg_stats_scan8 << "\t" << pct(time_sg_stats_ref, time_sg_stats_scan8) << endl;
#endif

    cout << "sw stats reference\t" << time_sw_stats_ref << endl;
#if ENABLE_WOZNIAK && HAVE_SSE2
    cout << "sw stats wozniak  \t" << time_sw_stats_woz << "\t" << pct(time_sw_stats_ref, time_sw_stats_woz) << endl;
#endif
#if ENABLE_STRIPED && HAVE_SSE2
    cout << "sw stats striped  \t" << time_sw_stats_striped << "\t" << pct(time_sw_stats_ref, time_sw_stats_striped) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE2
    cout << "sw stats scan     \t" << time_sw_stats_scan << "\t" << pct(time_sw_stats_ref, time_sw_stats_scan) << endl;
#endif
#if ENABLE_SCAN && HAVE_SSE41
    cout << "sw stats scan8    \t" << time_sw_stats_scan8 << "\t" << pct(time_sw_stats_ref, time_sw_stats_scan8) << endl;
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
