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
#include <tascel/UniformTaskCollIter.h>
#include <tascel/UniformTaskCollectionSplit.h>
#if THREADED
#include "tascelx.hpp"
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

#include "alignment.hpp"
#include "AlignStats.hpp"
#include "Bootstrap.hpp"
#include "combinations.h"
#include "constants.h"
#include "EdgeResult.hpp"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#if HAVE_ARMCI
#include "SequenceDatabaseArmci.hpp"
#endif
#include "SequenceDatabaseReplicated.hpp"
#include "Stats.hpp"
#define printf(...) fprintf(stdout, __VA_ARGS__); fflush(stdout);

using namespace std;
using namespace tascel;
using namespace pgraph;


typedef struct {
    unsigned long id;
} task_description_iter;

void createTaskFn(void *tsk, int tsk_size, TaskIndex tidx) {
    assert(tsk_size == sizeof(task_description_iter));
    *(unsigned long*)tsk = tidx;
}


typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_description;


typedef struct {
    AlignStats *stats;
    SequenceDatabase *sequences;
    vector<EdgeResult> *edge_results;
    Parameters *parameters;
    cell_t ***tbl;
    int ***del;
    int ***ins;
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


static void align(
    unsigned long seq_id[2],
    local_data_t *local_data,
    int thd)
{
    bool is_edge_answer = false;
    double t = 0;
    double tt = 0;
    int sscore;
    size_t max_len;

    AlignStats *stats = local_data->stats;
    SequenceDatabase *sequences = local_data->sequences;
    vector<EdgeResult> *edge_results = local_data->edge_results;
    Parameters *parameters = local_data->parameters;
    cell_t ***tbl = local_data->tbl;
    int ***del = local_data->del;
    int ***ins = local_data->ins;

    int open = parameters->open;
    int gap = parameters->gap;
    int AOL = parameters->AOL;
    int SIM = parameters->SIM;
    int OS = parameters->OS;

    tt = MPI_Wtime();

    Sequence &s1 = (*sequences)[seq_id[0]];
    Sequence &s2 = (*sequences)[seq_id[1]];
    unsigned long s1Len = s1.get_sequence_length();
    unsigned long s2Len = s2.get_sequence_length();
    bool do_alignment = true;
    if (parameters->use_length_filter) {
        int cutOff = AOL * SIM;
        if ((s1Len <= s2Len && (100 * s1Len < cutOff * s2Len))
                || (s2Len < s1Len && (100 * s2Len < cutOff * s1Len))) {
            stats[thd].work_skipped += s1Len * s2Len;
            ++stats[thd].align_skipped;
            do_alignment = false;
        }
    }

    if (do_alignment)
    {
        stats[thd].work += s1Len * s2Len;
        ++stats[thd].align_counts;
        t = MPI_Wtime();
        cell_t result = align_semi_affine(
                s1, s2, open, gap, tbl[thd], del[thd], ins[thd]);
        is_edge_answer = is_edge(
                result, s1, s2, AOL, SIM, OS, sscore, max_len);

        if (parameters->output_to_disk
                && (is_edge_answer || parameters->output_all))
        {
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
                        1.0*result.length/max_len,
                        1.0*result.matches/result.length,
                        1.0*result.score/sscore,
                        is_edge_answer));
        }
        if (is_edge_answer) {
            ++stats[thd].edge_counts;
        }
        t = MPI_Wtime() - t;
        stats[thd].time_align.push_back(t);
    }

    tt = MPI_Wtime() - tt;
    stats[thd].time_total += tt;
}


static void alignment_task_iter(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description_iter *desc = (task_description_iter*)bigd;
    local_data_t *local_data = (local_data_t*)pldata;
    unsigned long seq_id[2];
    AlignStats *stats = local_data->stats;
    double t = 0;

    t = MPI_Wtime();
    k_combination2(desc->id, seq_id);
    t = MPI_Wtime() - t;
    stats[thd].time_kcomb += t;

    align(seq_id, local_data, thd);
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
    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;

    align(seq_id, local_data, thd);
}


unsigned long populate_tasks_iter(
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

    --upper_limit; /* inclusive range */
    utcs[worker] = new UniformTaskCollIter(props, createTaskFn, lower_limit, upper_limit, worker);
    return upper_limit-lower_limit+1;
}


unsigned long populate_tasks(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        unsigned long ntasks, unsigned long tasks_per_worker, int worker,
        AlignStats *stats)
{
    int wrank = trank(worker);
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = lower_limit + tasks_per_worker;
    unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);
    double t;

    /* if I'm the last worker, add the remainder of the tasks */
    if (wrank == nprocs*NUM_WORKERS-1) {
        upper_limit += remainder;
    }

    task_description desc;
    unsigned long count = 0;
    unsigned long i;
    unsigned long seq_id[2];
    utcs[worker] = new UniformTaskCollectionSplit(props, worker);


    t = MPI_Wtime();
    k_combination(lower_limit, 2, seq_id);
    t = MPI_Wtime() - t;
    stats[worker].time_kcomb += t;
    for (i=lower_limit; i<upper_limit; ++i) {
        desc.id1 = seq_id[0];
        desc.id2 = seq_id[1];
        t = MPI_Wtime();
        next_combination(2, seq_id);
        t = MPI_Wtime() - t;
        stats[worker].time_kcomb += t;
        utcs[worker]->addTask(&desc, sizeof(desc));
        ++count;
    }

    return count;
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
    if (all_argv.size() <= 2 || all_argv.size() >= 4) {
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

    /* print parameters */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "-- Parameters ";
        cout << header.str() << endl;
        cout << *parameters << endl;
    }

    if (parameters->distribute_sequences) {
#if HAVE_ARMCI
        sequences = new SequenceDatabaseArmci(all_argv[1],
                parameters->memory_sequences, pgraph::comm, '\0');
#else
        if (0 == rank) {
            cerr << "config file specified to distribute sequences, "
                "but no implementation exists" << endl;
        }
        TascelConfig::finalize();
        pgraph::finalize();
        return 1;
#endif
    }
    else {
        sequences = new SequenceDatabaseReplicated(all_argv[1],
                parameters->memory_sequences, pgraph::comm, '\0');
    }

    local_data->sequences = sequences;
    local_data->tbl = new cell_t**[NUM_WORKERS];
    local_data->del = new int**[NUM_WORKERS];
    local_data->ins = new int**[NUM_WORKERS];
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        local_data->tbl[worker] = allocate_cell_table(2, sequences->longest());
        local_data->del[worker] = allocate_int_table(2, sequences->longest());
        local_data->ins[worker] = allocate_int_table(2, sequences->longest());
    }

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences->size(), 2);
    if (0 == rank) {
        cout << "brute force "
            << sequences->size()
            << " choose 2 has "
            << nCk
            << " combinations"
            << endl;
    }

    unsigned long ntasks = nCk;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 10;
    if (0 == rank) {
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(pgraph::comm);

    /* the tascel part */
    vector<double> populate_time(NUM_WORKERS, 0.0);
    vector<double> populate_count(NUM_WORKERS, 0.0);
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        TslFuncRegTbl *frt = NULL;
        TslFunc tf;
        TaskCollProps props;
        size_t taskSize = 0;

        frt = new TslFuncRegTbl;
        if (parameters->use_iterator) {
            tf = frt->add(alignment_task_iter);
            taskSize = sizeof(task_description_iter);
        }
        else {
            tf = frt->add(alignment_task);
            taskSize = sizeof(task_description);
        }

        props.functions(tf, frt)
            .taskSize(taskSize)
            .maxTasks(max_tasks_per_worker)
            .localData(local_data, sizeof(void*));
        populate_time[worker] = MPI_Wtime();
        if (parameters->use_iterator) {
            populate_count[worker] = populate_tasks_iter(
                    utcs, props, ntasks, tasks_per_worker, worker);
        }
        else {
            populate_count[worker] = populate_tasks(
                    utcs, props, ntasks, tasks_per_worker, worker, stats);
        }
        populate_time[worker] = MPI_Wtime() - populate_time[worker];
    }

    if (parameters->print_stats) {
        mpix_gather(populate_time, 0, pgraph::comm);
        mpix_gather(populate_count, 0, pgraph::comm);
        if (0 == rank) {
            Stats times;
            Stats counts;
            int p = cout.precision();
            char f = cout.fill();
            cout << setfill('-');
            cout << left;
            cout << setw(79) << "--- Task Population Statistics" << endl;
            cout << setfill(f);
            cout << setprecision(2);
            cout << right;
            cout << setw(5) << "pid";
            cout << setw(15) << "populate_time";
            cout << setw(15) << "populate_count";
            cout << endl;
            for (size_t i=0; i<populate_time.size(); ++i) {
                times.push_back(populate_time[i]);
                counts.push_back(populate_count[i]);
                cout << setw(5) << std::right << i
                    << setw(15) << fixed << populate_time[i] 
                    << setw(15) << fixed << populate_count[i]
                    << endl;
            }
            cout << string(79, '=') << endl;
            cout << "       " << Stats::header() << endl;
            cout << "Times  " << times << endl;
            cout << "Counts " << counts << endl;
            cout.precision(p);
        }
    }
    populate_time.clear();
    populate_count.clear();

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

    if (parameters->print_stats) {
        AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
        MPI_Gather(stats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
                rstats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
                0, pgraph::comm);
        /* synchronously print alignment stats all from process 0 */
        if (0 == rank) {
            Stats edge_counts;
            Stats align_counts;
            Stats align_skipped;
            Stats time_align;
            Stats time_kcomb;
            Stats time_total;
            Stats work;
            Stats work_skipped;
            int p = cout.precision();
            cout << setprecision(2);
            cout << right << setw(5) << "pid" << AlignStats::header() << endl;
            for(int i=0; i<nprocs*NUM_WORKERS; i++) {
                edge_counts.push_back(rstats[i].edge_counts);
                align_counts.push_back(rstats[i].align_counts);
                align_skipped.push_back(rstats[i].align_skipped);
                time_align.push_back(rstats[i].time_align.sum());
                time_kcomb.push_back(rstats[i].time_kcomb);
                time_total.push_back(rstats[i].time_total);
                work.push_back(rstats[i].work);
                work_skipped.push_back(rstats[i].work_skipped);
                cout << right << setw(5) << i << rstats[i] << endl;
            }
            cout << setprecision(1);
            cout << string(79, '=') << endl;
            cout << "           " << Stats::header() << endl;
            cout << "      Edges" << edge_counts << endl;
            cout << " Alignments" << align_counts << endl;
            cout << "  AlignSkip" << align_skipped << endl;
            cout << "     TAlign" << time_align << endl;
            cout << "     TTotal" << time_total << endl;
            cout << "       Work" << work << endl;
            cout << "WorkSkipped" << work_skipped << endl;
            cout.precision(p);
        }
        delete [] rstats;
    }

    if (parameters->print_stats) {
        StealingStats *stt = new StealingStats[NUM_WORKERS];
        StealingStats * rstt = new StealingStats[NUM_WORKERS*nprocs];

        for(int i=0; i<NUM_WORKERS; i++) {
            stt[i] = utcs[i]->getStats();
        }

        MPI_Gather(stt, sizeof(StealingStats)*NUM_WORKERS, MPI_CHAR, 
                rstt, sizeof(StealingStats)*NUM_WORKERS, MPI_CHAR, 
                0, pgraph::comm);

        /* synchronously print stealing stats all from process 0 */
        if (0 == rank) {
            StealingStats tstt;
            cout<<" pid"<<rstt[0].formatString()<<endl;      
            for(int i=0; i<nprocs*NUM_WORKERS; i++) {
                tstt += rstt[i];
                cout<<std::setw(4)<<std::right<<i<<rstt[i]<<endl;
            }
            cout << string(79, '=') << endl;
            cout<<"TOT "<<tstt<<endl;
        }

        delete [] rstt;
        delete [] stt;
    }

    amBarrier();
    MPI_Barrier(pgraph::comm);

    if (parameters->output_to_disk) {
        ofstream edge_out(get_edges_filename(rank).c_str());
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            for (size_t i=0,limit=edge_results[worker].size(); i<limit; ++i) {
                edge_out << edge_results[worker][i] << endl;
            }
        }
        edge_out.close();
    }

    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        delete utcs[worker];
        free_cell_table(local_data->tbl[worker], 2);
        free_int_table(local_data->del[worker], 2);
        free_int_table(local_data->ins[worker], 2);
    }
    delete [] utcs;
    delete [] stats;
    delete [] edge_results;
    delete parameters;
    delete sequences;
    delete [] local_data->tbl;
    delete [] local_data->del;
    delete [] local_data->ins;
    delete local_data;

    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}
