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
#include "EdgeResult.hpp"
#include "mpix.hpp"
#include "PairCheck.hpp"
#include "PairCheckGlobal.hpp"
#include "PairCheckLocal.hpp"
#include "PairCheckSemiLocal.hpp"
#include "PairCheckSmp.hpp"
#include "Parameters.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"
#include "Stats.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixBucketsTascel.hpp"
#include "SuffixTree.hpp"
#include "TreeStats.hpp"

#define ENABLE_ARMCI 0
#if HAVE_ARMCI && ENABLE_ARMCI
#include "SequenceDatabaseArmci.hpp"
#include "SuffixBucketsArmci.hpp"
#endif

#define printf(...) fprintf(stdout, __VA_ARGS__); fflush(stdout);

using namespace std;
using namespace tascel;
using namespace pgraph;


typedef struct {
    unsigned long id;
} task_description_one;

void createTaskFn(void *tsk, int tsk_size, TaskIndex tidx) {
    assert(tsk_size == sizeof(task_description_one));
    *(unsigned long*)tsk = tidx;
}


typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_description_two;


typedef struct {
    UniformTaskCollection **utcs;
    AlignStats *stats_align;
    TreeStats *stats_tree;
    SequenceDatabase *sequences;
    vector<EdgeResult> *edge_results;
    Parameters *parameters;
    PairCheck **pair_check;
    SuffixBuckets *suffix_buckets;
    cell_t ***tbl;
    int ***del;
    int ***ins;
} local_data_t;


#if HAVE_ARMCI && ENABLE_ARMCI
class DynamicTaskCounter {
    private:
        int rank_world;
        int size_world;
        unsigned long limit;
        local_data_t *local_data;
        void (*function)(unsigned long, local_data_t*, int);
        long **counter;
    public:
        DynamicTaskCounter(
                int rank_world, int size_world,
                unsigned long limit, local_data_t *local_data,
                void (*function)(unsigned long, local_data_t*, int))
            :   rank_world(rank_world)
            ,   size_world(size_world)
            ,   limit(limit)
            ,   local_data(local_data)
            ,   function(function)
            ,   counter(NULL)
        {
            initialize_armci();
        }
        void process() {
            long task_id = rank_world;
            long bytes = (0 == rank_world) ? sizeof(long) : 0;
            long **counter = new long*[size_world];
            int retval = ARMCI_Malloc((void**)counter, bytes);
            assert(0 == retval);
            if (0 == rank_world) {
                counter[0][0] = size_world; /* init counter */
            }
            ARMCI_Barrier();
            // TODO TAKS COUNTER AND ALIGNMENT LOOP
            while ((unsigned long)(task_id) < limit) {
                function((unsigned long)(task_id), local_data, 0);
                // next task
                retval = ARMCI_Rmw(ARMCI_FETCH_AND_ADD_LONG, &task_id, counter[0], 1, 0);
                assert(0 == retval);
            }
            retval = ARMCI_Free(counter[rank_world]);
            delete [] counter;
        }
        
};
#endif


static int trank(int thd);
static string get_edges_filename(int rank);
static void align(unsigned long seq_id[2], local_data_t *local_data, int thd);
static void alignment_task_iter_pre(
        unsigned long id, local_data_t *local_data, int thd);
static void alignment_task_iter(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd);
static void alignment_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd);
static unsigned long process_tree(unsigned long bid, Bucket *bucket, local_data_t *local_data, int thd);
static void tree_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd);
static void hybrid_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd);
static unsigned long populate_tasks_iter(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        unsigned long ntasks, unsigned long tasks_per_worker,
        int worker);
static unsigned long populate_tasks(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        unsigned long ntasks, unsigned long tasks_per_worker, AlignStats *stats,
        int worker);
static unsigned long populate_tasks_tree(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        local_data_t *local_data, int worker);
static unsigned long populate_tasks_tree_dynamic(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        local_data_t *local_data, int worker);
static unsigned long populate_tasks_tree_hybrid(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        local_data_t *local_data, int worker);


int main(int argc, char **argv)
{
    vector<string> all_argv;
    UniformTaskCollection **utcs = NULL;
    AlignStats *stats_align = NULL;
    TreeStats *stats_tree = NULL;
    SequenceDatabase *sequences = NULL;
    vector<EdgeResult> *edge_results = NULL;
    PairCheck **pair_check = NULL;
    SuffixBuckets *suffix_buckets = NULL;
    Parameters *parameters = NULL;
    local_data_t *local_data = NULL;
    char delimiter = '\0';

    pgraph::initialize(&argc, &argv);

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, pgraph::comm);

    /* initialize global data */
    utcs = new UniformTaskCollection*[NUM_WORKERS];
    stats_align = new AlignStats[NUM_WORKERS];
    stats_tree = new TreeStats[NUM_WORKERS];
    edge_results = new vector<EdgeResult>[NUM_WORKERS];
    pair_check = new PairCheck*[NUM_WORKERS];
    parameters = new Parameters;
    local_data = new local_data_t;
    local_data->utcs = utcs;
    local_data->stats_align = stats_align;
    local_data->stats_tree = stats_tree;
    local_data->edge_results = edge_results;
    local_data->pair_check = pair_check;
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

    if (parameters->use_tree
            || parameters->use_tree_dynamic
            || parameters->use_tree_hybrid) {
        delimiter = parameters->alphabet_dollar; /* dollar */
        if (parameters->dup_smp) {
            pair_check[0] = new PairCheckSmp;
            for (int worker=1; worker<NUM_WORKERS; ++worker) {
                pair_check[worker] = pair_check[0];
            }
        }
        else {
            for (int worker=0; worker<NUM_WORKERS; ++worker) {
                if (parameters->dup_local) {
                    pair_check[worker] = new PairCheckLocal;
                }
                else if (parameters->dup_semilocal) {
                    pair_check[worker] = new PairCheckSemiLocal;
                }
                else {
                    assert(0);
                }
            }
        }
    }

    if (parameters->distribute_sequences) {
#if HAVE_ARMCI && ENABLE_ARMCI
        sequences = new SequenceDatabaseArmci(all_argv[1],
                parameters->memory_sequences, pgraph::comm, delimiter);
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
                parameters->memory_sequences, pgraph::comm, delimiter);
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

    if (parameters->use_tree
            || parameters->use_tree_dynamic
            || parameters->use_tree_hybrid) {
#if HAVE_ARMCI && ENABLE_ARMCI
        suffix_buckets = new SuffixBucketsArmci(sequences, *parameters, pgraph::comm);
#else
        suffix_buckets = new SuffixBucketsTascel(sequences, *parameters, pgraph::comm);
#endif
        local_data->suffix_buckets = suffix_buckets;
    }

    /* how many combinations of sequences are there? */
    unsigned long ntasks = binomial_coefficient(sequences->size(), 2);
    if (0 == rank) {
        cout << "brute force "
            << sequences->size()
            << " choose 2 has "
            << ntasks
            << " combinations"
            << endl;
    }

    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 100;
    if (0 == rank) {
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(pgraph::comm);

    /* the tascel part */
    vector<double> populate_time(NUM_WORKERS, 0.0);
    vector<double> populate_count(NUM_WORKERS, 0.0);

    if (parameters->use_iterator) {
        TslFuncRegTbl *frt = new TslFuncRegTbl;
        TslFunc tf = frt->add(alignment_task_iter);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description_one))
            .maxTasks(parameters->memory_worker / sizeof(task_description_one))
            .localData(local_data, sizeof(void*));
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime();
            populate_count[worker] = populate_tasks_iter(
                    utcs, props, ntasks, tasks_per_worker, worker);
            populate_time[worker] = MPI_Wtime() - populate_time[worker];
        }
    }
    else if (parameters->use_tree) {
        TslFuncRegTbl *frt = new TslFuncRegTbl;
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description_two))
            .maxTasks(parameters->memory_worker / sizeof(task_description_two))
            .localData(local_data, sizeof(void*));
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime();
            populate_count[worker] = populate_tasks_tree(
                    utcs, props, local_data, worker);
            populate_time[worker] = MPI_Wtime() - populate_time[worker];
        }
    }
    else if (parameters->use_tree_dynamic) {
        TslFuncRegTbl *frt = new TslFuncRegTbl;
        TslFunc tf = frt->add(alignment_task);
        TslFunc tf2 = frt->add(tree_task);
        TaskCollProps props;
        props.functions(tf, frt, tf2)
            .taskSize(sizeof(task_description_two))
            .maxTasks(parameters->memory_worker / sizeof(task_description_two))
            .taskSize2(sizeof(task_description_one))
            .maxTasks2(suffix_buckets->owns().size())
            .localData(local_data, sizeof(void*));
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime();
            populate_count[worker] = populate_tasks_tree_dynamic(
                    utcs, props, local_data, worker);
            populate_time[worker] = MPI_Wtime() - populate_time[worker];
        }
    }
    else if (parameters->use_tree_hybrid) {
        TslFuncRegTbl *frt = new TslFuncRegTbl;
        TslFunc tf = frt->add(hybrid_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description_two))
            .maxTasks(parameters->memory_worker / sizeof(task_description_two))
            .localData(local_data, sizeof(void*));
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime();
            populate_count[worker] = populate_tasks_tree_hybrid(
                    utcs, props, local_data, worker);
            populate_time[worker] = MPI_Wtime() - populate_time[worker];
        }
    }
    else if (parameters->use_counter) {
    }
    else {
        TslFuncRegTbl *frt = new TslFuncRegTbl;
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description_two))
            .maxTasks(parameters->memory_worker / sizeof(task_description_two))
            .localData(local_data, sizeof(void*));
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime();
            populate_count[worker] = populate_tasks(
                    utcs, props, ntasks, tasks_per_worker, stats_align, worker);
            populate_time[worker] = MPI_Wtime() - populate_time[worker];
        }
    }

    if (suffix_buckets && parameters->print_stats) {
        if (0 == rank) {
            char f = cout.fill();
            cout << setfill('-');
            cout << left;
            cout << setw(79) << "--- SuffixBucket Statistics" << endl;
            cout << setfill(f);
            cout << "n_buckets " << suffix_buckets->size() << endl;
            cout << "non-empty " << suffix_buckets->size_nonempty() << endl;
        }
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

    if (parameters->use_counter) {
#if HAVE_ARMCI && ENABLE_ARMCI
        DynamicTaskCounter counter(
                rank, nprocs, ntasks, local_data, alignment_task_iter_pre);
        counter.process();
#else
        if (0 == rank) {
            cerr << "config file specified to use task counter, "
                "but no implementation exists" << endl;
        }
        TascelConfig::finalize();
        pgraph::finalize();
        return 1;
#endif
    }
    else {
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
    }

    amBarrier();
    MPI_Barrier(pgraph::comm);

    if (parameters->print_stats) {
        AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
        MPI_Gather(stats_align, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
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
            ostringstream header;
            int p = cout.precision();

            header.fill('-');
            header << left << setw(79) << "--- Align Stats ";
            cout << header.str() << endl;
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
            cout << string(79, '-') << endl;
            cout.precision(p);
        }
        delete [] rstats;
    }

    if (suffix_buckets && parameters->print_stats) {
        TreeStats * rstats = new TreeStats[NUM_WORKERS*nprocs];
        MPI_Gather(stats_tree, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
                rstats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
                0, pgraph::comm);
        /* synchronously print alignment stats all from process 0 */
        if (0 == rank) {
            TreeStats cumulative;
            ostringstream header;
            int p = cout.precision();

            header.fill('-');
            header << left << setw(79) << "--- Tree Stats ";
            cout << header.str() << endl;
            cout << setprecision(2);
            for(int i=0; i<nprocs*NUM_WORKERS; i++) {
                cout << right << setw(5) << i;
                cout << right << setw(14) << "name";
                cout << Stats::header() << endl;
                cumulative += rstats[i];
                cout << rstats[i] << endl;
            }
            cout << setprecision(1);
            cout << string(79, '=') << endl;
            cout << right << setw(5) << "TOTAL";
            cout << right << setw(14) << "name";
            cout << Stats::header() << endl;
            cout << cumulative << endl;
            cout << "first tree" << setw(25) << cumulative.time_first << endl;
            cout << " last tree" << setw(25) << cumulative.time_last << endl;
            cout << "      diff" << setw(25) << cumulative.time_last - cumulative.time_first << endl;
            cout.precision(p);
            cout << string(79, '-') << endl;
        }
        delete [] rstats;
    }

    if (!parameters->use_counter && parameters->print_stats) {
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
        if (!parameters->use_counter) {
            delete utcs[worker];
        }
        free_cell_table(local_data->tbl[worker], 2);
        free_int_table(local_data->del[worker], 2);
        free_int_table(local_data->ins[worker], 2);
    }
    delete [] utcs;
    delete [] stats_align;
    delete [] stats_tree;
    delete [] edge_results;
    if (parameters->use_tree
            || parameters->use_tree_dynamic
            || parameters->use_tree_hybrid) {
        int limit = parameters->dup_smp ? 1 : NUM_WORKERS;
        for (int worker=0; worker<limit; ++worker) {
            delete pair_check[worker];
        }
    }
    delete [] pair_check;
    delete suffix_buckets;
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

    AlignStats *stats = local_data->stats_align;
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
        do_alignment = SuffixTree::length_filter(s1Len, s2Len, AOL*SIM/100);
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
    else {
        stats[thd].work_skipped += s1Len * s2Len;
        stats[thd].align_skipped += 1;
    }

    tt = MPI_Wtime() - tt;
    stats[thd].time_total += tt;
}


static void alignment_task_iter_pre(
        unsigned long id, local_data_t *local_data, int thd)
{
    unsigned long seq_id[2];
    AlignStats *stats= local_data->stats_align;
    double t = 0;

    t = MPI_Wtime();
    k_combination2(id, seq_id);
    t = MPI_Wtime() - t;
    stats[thd].time_kcomb += t;

    align(seq_id, local_data, thd);
}


static void alignment_task_iter(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description_one *desc = (task_description_one*)bigd;
    local_data_t *local_data = (local_data_t*)pldata;
    alignment_task_iter_pre(desc->id, local_data, thd);
}


static void alignment_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description_two *desc = (task_description_two*)bigd;
    local_data_t *local_data = (local_data_t*)pldata;
    unsigned long seq_id[2];
    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;

    align(seq_id, local_data, thd);
}


static unsigned long process_tree(unsigned long bid, Bucket *bucket, local_data_t *local_data, int worker)
{
    unsigned long count = 0;
    SequenceDatabase *sequences = local_data->sequences;
    UniformTaskCollection **utcs = local_data->utcs;
    Parameters *parameters = local_data->parameters;
    TreeStats *stats = local_data->stats_tree;
    PairCheck **pair_check = local_data->pair_check;

    if (NULL != bucket->suffixes) {
        double t;
        set<pair<size_t,size_t> > local_pairs;

        assert(bucket->size > 0);

        if (stats[worker].time_first == 0.0) {
            stats[worker].time_first = MPI_Wtime();
        }

        /* construct tree */
        t = MPI_Wtime();
        SuffixTree *tree = new SuffixTree(sequences, bucket, *parameters);
        stats[worker].time_build.push_back(MPI_Wtime() - t);

        /* gather stats */
        stats[worker].trees++;
        stats[worker].suffixes.push_back(bucket->size);
        stats[worker].size.push_back(tree->get_size());
        stats[worker].size_internal.push_back(tree->get_size_internal());
        stats[worker].fanout.push_back(tree->get_fanout());
        stats[worker].depth.push_back(tree->get_depth());
        stats[worker].suffix_length.push_back(tree->get_suffix_length());

        /* generate pairs */
        t = MPI_Wtime();
        tree->generate_pairs(local_pairs);
        size_t orig_size = local_pairs.size();
        stats[worker].pairs.push_back(orig_size);
        stats[worker].time_process.push_back(MPI_Wtime() - t);
        delete tree;
        local_pairs = pair_check[worker]->check(local_pairs);
        size_t new_size = local_pairs.size();
#if 0
        cout << bid << "=" << local_data->suffix_buckets->bucket_kmer(bid)
            << "\t" << orig_size << "-" << new_size << "=" << (orig_size-new_size)
            << "\t" << utcs[worker]->size() << endl;
#endif

        /* populate task queue */
        for (set<pair<size_t,size_t> >::iterator it=local_pairs.begin();
                it!=local_pairs.end(); ++it) {
            task_description_two desc;
            desc.id1 = it->first;
            desc.id2 = it->second;
            utcs[worker]->addTask(&desc, sizeof(desc));
            ++count;
        }

        stats[worker].time_last = MPI_Wtime();
    }

    return count;
}


static void tree_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description_one *tree_desc = (task_description_one*)bigd;
    local_data_t *local_data = (local_data_t*)pldata;
    unsigned long i = tree_desc->id;
    Bucket *bucket = local_data->suffix_buckets->get(i);
    process_tree(i, bucket, local_data, thd);
    local_data->suffix_buckets->rem(bucket);
}


static void hybrid_task(
        UniformTaskCollection * utc,
        void *bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs,
        int thd)
{
    task_description_two *desc = (task_description_two*)bigd;
    if (desc->id1 == desc->id2) {
        tree_task(utc, bigd, bigd_len, pldata, pldata_len, data_bufs, thd);
    }
    else {
        alignment_task(utc, bigd, bigd_len, pldata, pldata_len, data_bufs, thd);
    }
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
        unsigned long ntasks, unsigned long tasks_per_worker, AlignStats *stats,
        int worker)
{
    int wrank = trank(worker);
    int wsize = nprocs*NUM_WORKERS;
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = (wrank+1)*tasks_per_worker;
    unsigned long remainder = ntasks % wsize;
    double t;

    if (wrank < remainder) {
        lower_limit += wrank;
        upper_limit += wrank+1;
    }
    else {
        lower_limit += remainder;
        upper_limit += remainder;
    }
    if (upper_limit > ntasks) {
        upper_limit = ntasks;
    }

    task_description_two desc;
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


static unsigned long populate_tasks_tree(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        local_data_t *local_data, int worker)
{
    SuffixBuckets *suffix_buckets = local_data->suffix_buckets;
    unsigned long count=0;
    const vector<size_t> &my_buckets = suffix_buckets->owns();

    utcs[worker] = new UniformTaskCollectionSplit(props, worker);

    for (size_t i=0; i<my_buckets.size(); ++i) {
        if ((i%size_t(NUM_WORKERS)) == size_t(worker)) {
            Bucket *bucket = suffix_buckets->get(my_buckets[i]);
            count += process_tree(my_buckets[i], bucket, local_data, worker);
            local_data->suffix_buckets->rem(bucket);
        }
    }

    return count;
}


static unsigned long populate_tasks_tree_dynamic(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        local_data_t *local_data, int worker)
{
    SuffixBuckets *suffix_buckets = local_data->suffix_buckets;
    unsigned long count=0;
    const vector<size_t> &my_buckets = suffix_buckets->owns();

    utcs[worker] = new UniformTaskCollectionSplit(props, worker);

    for (size_t i=0; i<my_buckets.size(); ++i) {
        if ((i%size_t(NUM_WORKERS)) == size_t(worker)) {
            size_t bid = my_buckets[i];
            Bucket *bucket = suffix_buckets->get(bid);
            if (NULL != bucket->suffixes) {
                task_description_one desc;
                desc.id = bid;
                utcs[worker]->addTask2(&desc, sizeof(desc));
                ++count;
            }
            local_data->suffix_buckets->rem(bucket);
        }
    }

    return count;
}


static unsigned long populate_tasks_tree_hybrid(
        UniformTaskCollection **utcs, const TaskCollProps &props,
        local_data_t *local_data, int worker)
{
    SuffixBuckets *suffix_buckets = local_data->suffix_buckets;
    unsigned long count=0;
    const vector<size_t> &my_buckets = suffix_buckets->owns();

    utcs[worker] = new UniformTaskCollectionSplit(props, worker);

    for (size_t i=0; i<my_buckets.size(); ++i) {
        if ((i%size_t(NUM_WORKERS)) == size_t(worker)) {
            size_t bid = my_buckets[i];
            Bucket *bucket = suffix_buckets->get(bid);
            if (NULL != bucket->suffixes) {
                task_description_two desc;
                desc.id1 = bid;
                desc.id2 = bid;
                utcs[worker]->addTask(&desc, sizeof(desc));
                ++count;
            }
            local_data->suffix_buckets->rem(bucket);
        }
    }

    return count;
}


