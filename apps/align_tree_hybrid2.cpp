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
#include <tascel/UniformTaskCollectionSplit.h>
#include <tascel/tmpi.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
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
#include "EdgeResult.hpp"
#include "Bootstrap.hpp"
#include "TreeStats.hpp"
#include "constants.h"
#include "combinations.h"
#include "alignment.hpp"
#include "mpix.hpp"
#include "tascelx.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseArmci.hpp"
#include "SequenceDatabaseReplicated.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixTree.hpp"

#define DEBUG 0

using namespace std;
using namespace tascel;
using namespace pgraph;


UniformTaskCollectionSplit** utcs = 0;
AlignStats *stats = 0;
TreeStats *treestats = 0;
SequenceDatabase *sequences = 0;
vector<EdgeResult> *edge_results = 0;
Parameters parameters;
SuffixBuckets2 *suffix_buckets = NULL;


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


typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_desc_hybrid;


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
    task_desc_hybrid *desc = (task_desc_hybrid*)_bigd;
    unsigned long seq_id[2];
    cell_t result;
    bool is_edge_answer = false;
    double t = 0;
    int sscore;
    size_t max_len;

    int open = parameters.open;
    int gap = parameters.gap;
    int AOL = parameters.AOL;
    int SIM = parameters.SIM;
    int OS = parameters.OS;

    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
    {
        t = MPI_Wtime();
        unsigned long s1Len = (*sequences)[seq_id[0]].get_sequence_length();
        unsigned long s2Len = (*sequences)[seq_id[1]].get_sequence_length();
        stats[thd].work += s1Len * s2Len;
        ++stats[thd].align_counts;
        sequences->align(seq_id[0], seq_id[1], result.score, result.matches, result.length, open, gap, thd);
        is_edge_answer = sequences->is_edge(
                seq_id[0],
                seq_id[1],
                result.score, result.matches, result.length,
                AOL, SIM, OS,
                sscore, max_len);

        if (is_edge_answer) {
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
                        1.0*result.length/max_len,
                        1.0*result.matches/result.length,
                        1.0*result.score/sscore
                        ));
            if (is_edge_answer) {
                ++stats[thd].edge_counts;
            }
        }
        t = MPI_Wtime() - t;
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
    }
}


static void tree_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int worker) {
    task_desc_hybrid *tree_desc = (task_desc_hybrid*)_bigd;
    unsigned long i = tree_desc->id1;
#if DEBUG
    cout << trank(worker) << " tree_task id1=" << i << endl;
#endif
    set<pair<size_t,size_t> > local_pairs;
    /* suffix tree construction & processing */
    Bucket *bucket = suffix_buckets->get(i);
    if (NULL != bucket->suffixes) {
        double atimer = MPI_Wtime();
        SuffixTree *tree = new SuffixTree(
                sequences, bucket, parameters);
        treestats[worker].size_internal += tree->get_size_internal();
        treestats[worker].fanout += tree->get_avg_fanout();
        treestats[worker].avgdepth += tree->get_avg_depth();
        treestats[worker].deepest += tree->get_deepest();
        treestats[worker].suffix_avg_length += tree->get_suffix_avg_length();
        treestats[worker].suffix_max_length += tree->get_suffix_max_length();
        treestats[worker].time_build += MPI_Wtime() - atimer;
        treestats[worker].times += MPI_Wtime() - atimer;
        treestats[worker].trees += 1;
        treestats[worker].suffixes += bucket->size;
        atimer = MPI_Wtime();
        tree->generate_pairs(local_pairs);
        treestats[worker].time_process += MPI_Wtime() - atimer;
        treestats[worker].times += MPI_Wtime() - atimer;
        delete tree;
        suffix_buckets->rem(i, bucket);
    }
    
    task_desc_hybrid desc;
    for (set<pair<size_t,size_t> >::iterator it=local_pairs.begin();
            it!=local_pairs.end(); ++it) {
        desc.id1 = it->first;
        desc.id2 = it->second;
#if DEBUG
        cout << trank(worker) << " added " << desc.id1 << "," << desc.id2 << endl;
#endif
        utcs[worker]->addTask(&desc, sizeof(desc));
    }
}


static void hybrid_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_desc_hybrid *desc = (task_desc_hybrid*)_bigd;
    if (desc->id1 == desc->id2) {
        tree_task(utc, _bigd, bigd_len, pldata, pldata_len, data_bufs, thd);
    }
    else {
        alignment_task(utc, _bigd, bigd_len, pldata, pldata_len, data_bufs, thd);
    }
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    vector<string> all_argv;
    unsigned long nCk;
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */

    pgraph::initialize(&argc, &argv);

    /* initialize tascel */
    double totaltime = MPI_Wtime();
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, comm);
    utcs = new UniformTaskCollectionSplit*[NUM_WORKERS];
    stats = new AlignStats[NUM_WORKERS];
    treestats = new TreeStats[NUM_WORKERS];
    edge_results = new vector<EdgeResult>[NUM_WORKERS];

    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank) {
        MPI_CHECK(MPI_Bcast(&argc, 1, MPI_INT, 0, comm));
        for (int i=0; i<argc; ++i) {
            int length = strlen(argv[i])+1;
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(argv[i], length, MPI_CHAR, 0, comm));
            all_argv.push_back(argv[i]);
        }
    } else {
        int all_argc;
        MPI_CHECK(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm));
        for (int i=0; i<all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm));
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
        MPI_Barrier(comm);
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
            printf("usage: align sequence_file memory_budget param\n");
        }
        pgraph::finalize();
        return 1;
    }
    else if (all_argv.size() >= 4) {
        parameters.parse(all_argv[3].c_str(), comm);
    }

#define GB (1073741824U)
#if 1
    sequences = new SequenceDatabaseArmci(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, DOLLAR);
#else
    sequences = new SequenceDatabaseReplicated(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, DOLLAR);
#endif

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences->get_global_count(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences->get_global_count(), nCk);
    }

    (void) time(&t1);
    suffix_buckets = new SuffixBuckets2(sequences, parameters, comm);
    mpix_print_sync("after SuffixBuckets2 ctor", comm);
    (void) time(&t2);
    if (0 == trank(0)) {
        printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
        fflush(stdout);
    }

    unsigned long ntasks = nCk;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long max_tasks_per_worker = 2*(ntasks/global_num_workers + 1);
    max_tasks_per_worker += suffix_buckets->n_buckets/global_num_workers + 1;
    unsigned long GB_4 = 268435456U;
    unsigned long GB_2 = 536870912U;
    //max_tasks_per_worker = std::min(max_tasks_per_worker, GB_4/sizeof(task_desc_hybrid));
    max_tasks_per_worker = GB/sizeof(task_desc_hybrid);
    if (0 == trank(0)) {
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
        printf("                GB_2=%lu\n", GB_2);
        fflush(stdout);
    }
    if (0 == trank(0)) {
        printf("----------------------------------------------\n");
        printf("%-20s: %d\n", "slide size", parameters.window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters.exact_match_len);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters.AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters.SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters.OS);
        printf("----------------------------------------------\n");
        fflush(stdout);
    }

    MPI_Barrier(comm);

    /* the tascel part */
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        UniformTaskCollectionSplit*& utc = utcs[worker];
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(hybrid_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_desc_hybrid))
            .maxTasks(max_tasks_per_worker);
        utc = new UniformTaskCollectionSplit(props, worker);
    }
    if (0 == trank(0)) {
        printf("created UTCs\n");
        fflush(stdout);
    }

    int worker = 0;
    for (size_t j=rank; j<suffix_buckets->n_buckets; j+=nprocs) {
        int bucket_index = j / nprocs;
        if (suffix_buckets->buckets[bucket_index].size > 0) {
            task_desc_hybrid desc;
            desc.id1 = j;
            desc.id2 = j;
#if DEBUG
            cout << trank(worker) << " added " << desc.id1 << endl;
#endif
            utcs[worker]->addTask(&desc, sizeof(task_desc_hybrid));
            worker = (worker + 1) % NUM_WORKERS;
        }
    }
#if DEBUG || 1
    mpix_print_sync("done adding tasks", rank, comm);
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
    MPI_Barrier(comm);

#if 1
    size_t *gremote_buckets = new size_t[nprocs];
    MPI_Gather(&suffix_buckets->count_remote_buckets, sizeof(size_t), MPI_CHAR,
            gremote_buckets, sizeof(size_t), MPI_CHAR, 0, comm);
    if (0 == rank) {
        cout << setw(4) << right << "pid " << "remote-buckets" << endl;
        for(int i=0; i<nprocs; i++) {
            cout << std::setw(4) << right << i << " " << gremote_buckets[i] << endl;
        }
    }
    delete [] gremote_buckets;
    size_t *gremote_suffixes = new size_t[nprocs];
    MPI_Gather(&suffix_buckets->count_remote_suffixes, sizeof(size_t), MPI_CHAR,
            gremote_suffixes, sizeof(size_t), MPI_CHAR, 0, comm);
    if (0 == rank) {
        cout << setw(4) << right << "pid " << "remote-suffixes" << endl;
        for(int i=0; i<nprocs; i++) {
            cout << std::setw(4) << right << i << " " << gremote_suffixes[i] << endl;
        }
    }
    delete [] gremote_suffixes;
#endif

#if 1
    TreeStats * gtreestats = new TreeStats[NUM_WORKERS*nprocs];
    MPI_Gather(treestats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
	       gtreestats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print tree stats all from process 0 */
    if (0 == rank) {
        TreeStats totals;
        TreeStats mins = gtreestats[0];
        TreeStats maxs = gtreestats[0];
        cout << setw(4) << right << "pid " << gtreestats[0].getHeader() << endl;
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += gtreestats[i];
            mins < gtreestats[i];
            maxs > gtreestats[i];
            cout << std::setw(4) << right << i << " " << gtreestats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT " << totals << endl;
        cout << setw(4) << right << "MIN " << mins << endl;
        cout << setw(4) << right << "MAX " << maxs << endl;
        cout << setw(4) << right << "AVG " << (totals/(nprocs*NUM_WORKERS)) << endl;
        cout << setw(4) << right << "STD " << (totals.stddev(nprocs*NUM_WORKERS, gtreestats)) << endl;
        cout << "==============================================" << endl;
    }
    delete [] gtreestats;
    gtreestats=NULL;
#endif

    AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
    MPI_Gather(stats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       rstats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank) {
        AlignStats totals;
        AlignStats mins;
        AlignStats maxs;
        cout << setw(4) << right << "pid" << rstats[0].getHeader() << endl;
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += rstats[i];
            mins < rstats[i];
            maxs > rstats[i];
            cout << std::setw(4) << right << i << " " << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT " << totals << endl;
        cout << setw(4) << right << "MIN " << mins << endl;
        cout << setw(4) << right << "MAX " << maxs << endl;
        cout << setw(4) << right << "AVG " << (totals/(nprocs*NUM_WORKERS)) << endl;
        cout << setw(4) << right << "STD " << (totals.stddev(nprocs*NUM_WORKERS, rstats)) << endl;
        cout << "==============================================" << endl;
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
	       0, comm);

    /* synchronously print stealing stats all from process 0 */
    if (0 == rank) {
      StealingStats tstt;
      cout<<" pid "<<rstt[0].formatString()<<endl;      
      for(int i=0; i<nprocs*NUM_WORKERS; i++) {
	tstt += rstt[i];
	cout<<std::setw(4)<<std::right<<i<<" "<<rstt[i]<<endl;
      }
      cout<<"=============================================="<<endl;
      cout<<"TOT "<<tstt<<endl;
    }
    delete [] rstt;
    rstt=NULL;
    delete [] stt;

    amBarrier();
    MPI_Barrier(comm);

    ofstream edge_out(get_edges_filename(rank).c_str());
    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        delete utcs[worker];
        for (size_t i=0,limit=edge_results[worker].size(); i<limit; ++i) {
            edge_out << edge_results[worker][i] << endl;
        }
    }
    edge_out.close();

    delete suffix_buckets;
    delete [] utcs;
    delete [] edge_results;
    delete sequences;

    totaltime = MPI_Wtime() - totaltime;
    if (0 == trank(0)) {
        cout << "totaltime = " << totaltime << endl;
    }
    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}
