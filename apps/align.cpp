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
UniformTaskCollIter** utcs = 0;
#else
UniformTaskCollectionSplit** utcs = 0;
#endif
AlignStats *stats = 0;
SequenceDatabase *sequences = 0;
vector<EdgeResult> *edge_results = 0;
Parameters parameters;

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





static size_t parse_memory_budget(const string& value)
{
    long budget = 0;
    char budget_multiplier = 0;
    istringstream iss(value);

    if (isdigit(*value.rbegin())) {
        iss >> budget;
    }
    else {
        iss >> budget >> budget_multiplier;
    }

    if (budget <= 0) {
        cerr << "memory budget must be positive real number" << endl;
        assert(budget > 0);
    }

    if (budget_multiplier == 'b' || budget_multiplier == 'B') {
        budget *= 1; /* byte */
    }
    else if (budget_multiplier == 'k' || budget_multiplier == 'K') {
        budget *= 1024; /* kilobyte */
    }
    else if (budget_multiplier == 'm' || budget_multiplier == 'M') {
        budget *= 1048576; /* megabyte */
    }
    else if (budget_multiplier == 'g' || budget_multiplier == 'G') {
        budget *= 1073741824; /* gigabyte */
    }
    else if (budget_multiplier != 0) {
        cerr << "unrecognized size multiplier" << endl;
        assert(0);
    }

    assert(budget > 0);
    return size_t(budget);
}


static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_description *desc = (task_description*)_bigd;
    unsigned long seq_id[2];
    cell_t result;
    bool is_edge_answer = false;
    double t = 0;
    double tt = 0;
    int sscore;
    size_t max_len;

    int open = -10;
    int gap = -1;
    int AOL = 8;
    int SIM = 4;
    int OS = 3;

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


#ifdef USE_ITER
#else
unsigned long populate_tasks(
        unsigned long ntasks, unsigned long tasks_per_worker, int worker)
{
    task_description desc;
    int wrank = trank(worker);
    unsigned long count = 0;
    unsigned long i;
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = lower_limit + tasks_per_worker;
    unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);

    unsigned long seq_id[2];
    k_combination(lower_limit, 2, seq_id);
    for (i=lower_limit; i<upper_limit; ++i) {
        desc.id1 = seq_id[0];
        desc.id2 = seq_id[1];
#if DEBUG
        cout << wrank << " added " << seq_id[0] << "," << seq_id[1] << endl;
#endif
        next_combination(2, seq_id);
        utcs[worker]->addTask(&desc, sizeof(desc));
        count++;
    }
    /* if I'm the last worker, add the remainder of the tasks */
    if (wrank == nprocs*NUM_WORKERS-1) {
        for (/*ignore*/; i<upper_limit+remainder; ++i) {
            count++;
            desc.id1 = seq_id[0];
            desc.id2 = seq_id[1];
            next_combination(2, seq_id);
            utcs[worker]->addTask(&desc, sizeof(desc));
        }
    }

    return count;
}
#endif


int main(int argc, char **argv)
{
    vector<string> all_argv;
    unsigned long nCk;

    pgraph::initialize(&argc, &argv);

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, pgraph::comm);
#ifdef USE_ITER
    utcs = new UniformTaskCollIter*[NUM_WORKERS];
#else
    utcs = new UniformTaskCollectionSplit*[NUM_WORKERS];
#endif
    stats = new AlignStats[NUM_WORKERS];
    edge_results = new vector<EdgeResult>[NUM_WORKERS];

    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank) {
        MPI_CHECK(MPI_Bcast(&argc, 1, MPI_INT, 0, pgraph::comm));
        for (int i=0; i<argc; ++i) {
            int length = strlen(argv[i])+1;
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, pgraph::comm));
            MPI_CHECK(MPI_Bcast(argv[i], length, MPI_CHAR, 0, pgraph::comm));
            all_argv.push_back(argv[i]);
        }
    } else {
        int all_argc;
        MPI_CHECK(MPI_Bcast(&all_argc, 1, MPI_INT, 0, pgraph::comm));
        for (int i=0; i<all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, pgraph::comm));
            MPI_CHECK(MPI_Bcast(buffer, length, MPI_CHAR, 0, pgraph::comm));
            all_argv.push_back(buffer);
        }
    }

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
    if (all_argv.size() <= 2 || all_argv.size() >= 5) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() <= 2) {
                printf("missing memory budget\n");
            }
            else if (all_argv.size() >= 5) {
                printf("too many arguments\n");
            }
            printf("usage: align sequence_file memory_budget\n");
        }
        pgraph::finalize();
        return 1;
    }
    else if (all_argv.size() >= 4) {
        parameters.parse(all_argv[3].c_str(), pgraph::comm);
    }
    if (0 == trank(0)) {
        printf("----------------------------------------------\n");
        printf("%-20s: %d\n", "slide size", parameters.window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters.exact_match_len);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters.AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters.SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters.OS);
        printf("----------------------------------------------\n");
    }

#ifdef USE_GARRAY
    sequences = new SequenceDatabaseGArray(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), '\0');
#elif HAVE_ARMCI
    sequences = new SequenceDatabaseArmci(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), pgraph::comm, NUM_WORKERS, '\0');
#else
    sequences = new SequenceDatabaseReplicated(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), pgraph::comm, NUM_WORKERS, '\0');
#endif

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences->get_global_count(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences->get_global_count(), nCk);
    }

    double selectivity = 1.0;
    if (all_argv.size() == 3) {
      selectivity  = fabs(min(1.0,atof(all_argv[2].c_str())));
    }
    unsigned long nalignments = (long)(0.5+selectivity*nCk);
    unsigned long ntasks = nalignments;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 10;
    if (0 == trank(0)) {
        printf("selectivity=%lf\n", selectivity);
        printf("nalignments=%lu\n", nalignments);
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(pgraph::comm);

    /* the tascel part */
#ifdef USE_ITER
#else
    double populate_times[NUM_WORKERS];
    unsigned long count[NUM_WORKERS];
#endif
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        edge_results[worker].reserve(tasks_per_worker);
#ifdef USE_ITER
        UniformTaskCollIter*& utc = utcs[worker];
#else
        UniformTaskCollectionSplit*& utc = utcs[worker];
#endif
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description))
            .maxTasks(max_tasks_per_worker);
#ifdef USE_ITER
        int wrank = trank(worker);
        unsigned long lower_limit = wrank*tasks_per_worker;
        unsigned long upper_limit = lower_limit + tasks_per_worker;
        unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);
        /* if I'm the last worker, add the remainder of the tasks */
        if (wrank == nprocs*NUM_WORKERS-1) {
            upper_limit += remainder;
        }
        --upper_limit;
        utc = new UniformTaskCollIter(props, createTaskFn, lower_limit, upper_limit, worker);
        //if (trank(worker) == 0) {
        //    utc = new UniformTaskCollIter(props, createTaskFn, 0, ntasks-1, worker);
        //}
        //else {
        //    utc = new UniformTaskCollIter(props, createTaskFn, 0, -1, worker);
        //}
#else
        utc = new UniformTaskCollectionSplit(props, worker);
        /* add some tasks */
        populate_times[worker] = MPI_Wtime();
        count[worker] = populate_tasks(ntasks, tasks_per_worker, worker);
        populate_times[worker] = MPI_Wtime() - populate_times[worker];
#endif
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

    allocate_threads();
    initialize_threads(utcs);

    double mytimer = MPI_Wtime();
    utcs[0]->process();
    mytimer = MPI_Wtime() - mytimer;
    if (0 == trank(0)) {
        cout << "mytimer=" << mytimer << endl;
    }

    finalize_threads();

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
    delete [] stats;

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
    delete [] edge_results;
    delete sequences;

    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}
