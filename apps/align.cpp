/**
 * @author jeff.daily@pnnl.gov
 *
 * Alignment of an input dataset using work stealing. Each MPI task reads the
 * input file.
 */
#include "config.h"

#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollection.h>
#ifdef USE_ITER
#include <tascel/UniformTaskCollIter.h>
#else
#include <tascel/UniformTaskCollectionSplit.h>
#endif

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AlignStats.hpp"
#include "Bootstrap.hpp"
#include "constants.h"
#include "combinations.h"
#include "EdgeResult.hpp"
#include "alignment.hpp"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "tascelx.hpp"
#include "SequenceDatabase.hpp"
#ifdef USE_GARRAY
#include "SequenceDatabaseGArray.hpp"
#elif HAVE_ARMCI
#include "SequenceDatabaseArmci.hpp"
#else
#include "SequenceDatabaseReplicated.hpp"
#endif
#define printf(...) fprintf(stdout, __VA_ARGS__); fflush(stdout);

using namespace std;
using namespace tascel;
using namespace pgraph;

#define ALL_RESULTS 1

#ifdef USE_ITER
typedef struct {
    unsigned long id;
} task_description;
void createTaskFn(void *tsk, int tsk_size, TaskIndex tidx) {
    assert(tsk_size == sizeof(task_description));
    *(unsigned long*)tsk = tidx;
}
#else
typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_description;
#endif


typedef struct {
    AlignStats *stats;
    SequenceDatabase *sequences;
    vector<EdgeResult> *edge_results;
    Parameters *parameters;
} local_data_t;


static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}


static string get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}


static void alignment_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description *desc = (task_description*)bigd;
    local_data_t *local_data = (local_data_t*)pldata;
    unsigned long seq_id[2];
    cell_t result;
    bool is_edge_answer = false;
    double t = 0;
    double tt = 0;
    int sscore;
    size_t max_len;

    AlignStats *stats = local_data->stats;
    SequenceDatabase *sequences = local_data->sequences;
    vector<EdgeResult> *edge_results = local_data->edge_results;
    Parameters *parameters = local_data->parameters;

    int open = parameters->open;
    int gap = parameters->gap;
    int AOL = parameters->AOL;
    int SIM = parameters->SIM;
    int OS = parameters->OS;

    tt = MPI_Wtime();

#ifdef USE_ITER
    t = MPI_Wtime();
    k_combination2(desc->id, seq_id);
    t = MPI_Wtime() - t;
    stats[thd].kcomb_times_tot += t;
#else
    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
#endif
    unsigned long s1Len = (*sequences)[seq_id[0]].get_sequence_length();
    unsigned long s2Len = (*sequences)[seq_id[1]].get_sequence_length();
#ifdef LENGTH_FILTER
    int cutOff = AOL * SIM;
    if ((s1Len <= s2Len && (100 * s1Len < cutOff * s2Len))
            || (s2Len < s1Len && (100 * s2Len < cutOff * s1Len))) {
        stats[thd].work_skipped += s1Len * s2Len;
        ++stats[thd].align_skipped;
    }
    else
#endif
    {
        stats[thd].work += s1Len * s2Len;
        ++stats[thd].align_counts;
        t = MPI_Wtime();
#if defined(NOALIGN)
#else
#if USE_SSW
        sequences->align_ssw(seq_id[0], seq_id[1], result.score, result.matches, result.length, open, gap, thd);
#else
        sequences->align(seq_id[0], seq_id[1], result.score, result.matches, result.length, open, gap, thd);
#endif
        is_edge_answer = sequences->is_edge(
                seq_id[0],
                seq_id[1],
                result.score, result.matches, result.length,
                AOL, SIM, OS,
                sscore, max_len);

        if (is_edge_answer || ALL_RESULTS)
        {
#if DEBUG
            cout << trank(thd)
                << ": aligned " << seq_id[0] << " " << seq_id[1]
                << ": (score,ndig,alen)=("
                << result.score << ","
                << result.matches << ","
                << result.length << ")"
                << ": edge? " << is_edge_answer << endl;
#endif
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
#if 0
                        1.0*result.alen/max_len,
                        1.0*result.ndig/result.alen,
                        1.0*result.score/sscore
#else
                        result.length,
                        result.matches,
                        result.score
#endif
#if ALL_RESULTS
                        ,is_edge_answer
#endif
                        ));
            if (is_edge_answer) {
                ++stats[thd].edge_counts;
            }
        }
#endif
        t = MPI_Wtime() - t;
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
    }

    tt = MPI_Wtime() - tt;
    stats[thd].total_times_tot += tt;
}


unsigned long populate_tasks(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        unsigned long ntasks, unsigned long tasks_per_worker, int worker)
{
    int wrank = trank(worker);
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = lower_limit + tasks_per_worker;
    unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);

    /* if I'm the last worker, add the remainder of the tasks */
    if (wrank == nprocs*NUM_WORKERS-1) {
        upper_limit += remainder;
    }

#ifdef USE_ITER
    --upper_limit;
    utcs[worker] = new UniformTaskCollIter(props, createTaskFn, lower_limit, upper_limit, worker);
    return upper_limit-lower_limit;
#else
    task_description desc;
    unsigned long count = 0;
    unsigned long i;
    unsigned long seq_id[2];
    utcs[worker] = new UniformTaskCollectionSplit(props, worker);

    k_combination(lower_limit, 2, seq_id);
    for (i=lower_limit; i<upper_limit; ++i) {
        desc.id1 = seq_id[0];
        desc.id2 = seq_id[1];
#if DEBUG
        cout << wrank << " added " << seq_id[0] << "," << seq_id[1] << endl;
#endif
        next_combination(2, seq_id);
        utcs[worker]->addTask(&desc, sizeof(desc));
        ++count;
    }

    return count;
#endif
}


int main(int argc, char **argv)
{
    vector<string> all_argv;
    unsigned long nCk = 0;
    UniformTaskCollection **utcs = NULL;
    AlignStats *stats = NULL;
    SequenceDatabase *sequences = NULL;
    vector<EdgeResult> *edge_results = NULL;
    Parameters *parameters = NULL;
    local_data_t *local_data = NULL;

    pgraph::initialize(&argc, &argv);

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, pgraph::comm);

    /* initialize global data */
    utcs = new UniformTaskCollection*[NUM_WORKERS];
    stats = new AlignStats[NUM_WORKERS];
    edge_results = new vector<EdgeResult>[NUM_WORKERS];
    parameters = new Parameters;
    local_data = new local_data_t;
    local_data->stats = stats;
    local_data->edge_results = edge_results;
    local_data->parameters = parameters;


    /* MPI standard does not guarantee all procs receive argc and argv */
    mpix_bcast_argv(argc, argv, all_argv, pgraph::comm);

#if DEBUG
    /* print the command line arguments */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            for (size_t j=0; j<all_argv.size(); ++j) {
                printf("[%d] argv[%zd]=%s\n", rank, j, all_argv[j].c_str());
            }
        }
        MPI_Barrier(pgraph::comm);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 2 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() >= 4) {
                printf("too many arguments\n");
            }
            printf("usage: align sequence_file <config_file>\n");
        }
        pgraph::finalize();
        return 1;
    }
    else if (all_argv.size() >= 3) {
        parameters->parse(all_argv[2].c_str(), pgraph::comm);
    }
    if (0 == trank(0)) {
        printf("----------------------------------------------\n");
        printf("%-20s: %d\n", "slide size", parameters->window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters->exact_match_len);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters->AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters->SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters->OS);
        printf("%-20s: %d\n", "gap open", parameters->open);
        printf("%-20s: %d\n", "gap extend", parameters->gap);
        printf("%-20s: %zu\n", "mem worker", parameters->mem_worker);
        printf("%-20s: %zu\n", "mem sequences", parameters->mem_sequences);
        printf("----------------------------------------------\n");
    }

#ifdef USE_GARRAY
    sequences = new SequenceDatabaseGArray(all_argv[1],
            parameters->mem_sequences, '\0');
#elif HAVE_ARMCI
    sequences = new SequenceDatabaseArmci(all_argv[1],
            parameters->mem_sequences, pgraph::comm, NUM_WORKERS, '\0');
#else
    sequences = new SequenceDatabaseReplicated(all_argv[1],
            parameters->mem_sequences, pgraph::comm, NUM_WORKERS, '\0');
#endif

    local_data->sequences = sequences;

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences->get_global_count(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences->get_global_count(), nCk);
    }

    unsigned long nalignments = nCk;
    unsigned long ntasks = nalignments;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 10;
    if (0 == trank(0)) {
        printf("nalignments=%lu\n", nalignments);
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(pgraph::comm);

    /* the tascel part */
    double populate_times[NUM_WORKERS];
    unsigned long count[NUM_WORKERS];
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;

        props.functions(tf, frt)
            .taskSize(sizeof(task_description))
            .maxTasks(max_tasks_per_worker)
            .localData(local_data, sizeof(void*));
        populate_times[worker] = MPI_Wtime();
        count[worker] = populate_tasks(utcs, props, ntasks, tasks_per_worker, worker);
        populate_times[worker] = MPI_Wtime() - populate_times[worker];
    }
#if DEBUG
    double *g_populate_times = new double[nprocs*NUM_WORKERS];
    MPI_CHECK(MPI_Gather(populate_times, NUM_WORKERS, MPI_DOUBLE,
                g_populate_times, NUM_WORKERS, MPI_DOUBLE, 0, pgraph::comm));
    if (0 == rank) {
        double tally = 0;
        cout << " pid populate_time" << endl;
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            tally += g_populate_times[i];
            cout << std::setw(4) << std::right << i
                << setw(14) << fixed << g_populate_times[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "T" << setw(14) << fixed << tally << endl;
    }
    delete [] g_populate_times;
#endif
#if DEBUG
    unsigned long *task_counts = new unsigned long[nprocs*NUM_WORKERS];
    MPI_CHECK(MPI_Gather(count, NUM_WORKERS, MPI_UNSIGNED_LONG,
                task_counts, NUM_WORKERS, MPI_UNSIGNED_LONG, 0, pgraph::comm));
    if (0 == rank) {
        unsigned long tally = 0;
        cout << " pid        ntasks" << endl;
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            tally += task_counts[i];
            cout << std::setw(4) << std::right << i
                << setw(14) << task_counts[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "T" << setw(14) << tally << endl;
        if (tally != ntasks) {
            cout << "tally != ntasks\t" << tally << " != " << ntasks << endl;
        }
    }
    delete [] task_counts;
#endif

    amBarrier();

#if THREADED
    allocate_threads();
    initialize_threads(utcs);
#endif

    double mytimer = MPI_Wtime();
    utcs[0]->process();
    mytimer = MPI_Wtime() - mytimer;
    if (0 == trank(0)) {
        cout << "mytimer=" << mytimer << endl;
    }

#if THREADED
    finalize_threads();
#endif

    amBarrier();
    MPI_Barrier(pgraph::comm);

    AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
    MPI_Gather(stats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       rstats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       0, pgraph::comm);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank) {
        AlignStats totals;
        cout << " pid" << rstats[0].getHeader() << endl;      
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += rstats[i];
            cout << std::setw(4) << std::right << i << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT" << totals << endl;
    }
    delete [] rstats;
    rstats=NULL;

    StealingStats *stt = new StealingStats[NUM_WORKERS];
    for(int i=0; i<NUM_WORKERS; i++) {
        stt[i] = utcs[i]->getStats();
    }

    StealingStats * rstt = new StealingStats[NUM_WORKERS*nprocs];
    MPI_Gather(stt, sizeof(StealingStats)*NUM_WORKERS, MPI_CHAR, 
            rstt, sizeof(StealingStats)*NUM_WORKERS, MPI_CHAR, 
            0, pgraph::comm);

    /* synchronously print stealing stats all from process 0 */
    if (0 == rank) {
        StealingStats tstt;
        cout<<" pid "<<rstt[0].formatString()<<endl;      
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            tstt += rstt[i];
            cout<<std::setw(4)<<std::right<<i<<rstt[i]<<endl;
        }
        cout<<"=============================================="<<endl;
        cout<<"TOT "<<tstt<<endl;
    }
    delete [] rstt;
    rstt=NULL;
    delete [] stt;

    amBarrier();
    MPI_Barrier(pgraph::comm);

#if OUTPUT_EDGES
    ofstream edge_out(get_edges_filename(rank).c_str());
#endif
    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        delete utcs[worker];
#if OUTPUT_EDGES
        for (size_t i=0,limit=edge_results[worker].size(); i<limit; ++i) {
            edge_out << edge_results[worker][i] << endl;
        }
#endif
    }
#if OUTPUT_EDGES
    edge_out.close();
#endif

    delete [] utcs;
    delete [] stats;
    delete [] edge_results;
    delete parameters;
    delete sequences;
    delete local_data;

    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}
