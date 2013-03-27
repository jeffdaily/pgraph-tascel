/**
 * @author jeff.daily@pnnl.gov
 *
 * A first attempt at a brute force alignment of an input dataset using work
 * stealing. Each MPI task reads the input file.
 */
#include "config.h"

#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollection.h>
#include <tascel/UniformTaskCollIter.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <strstream>
#include <vector>

#include "AlignStats.h"
#include "combinations.h"
#include "dynamic.h"
#include "mpix.h"

using namespace std;
using namespace tascel;
using namespace pgraph;


int rank = 0;
int nprocs = 0;
cell_t ***tbl = 0;
int ***del = 0;
int ***ins = 0;
UniformTaskCollIter** utcs = 0;
AlignStats *stats = 0;
static pthread_t *threadHandles = 0;
static unsigned *threadRanks = 0;
vector<string> sequences;
// Synchronization for worker threads
pthread_barrier_t workersStart, workersEnd;
// Synchronization for server thread
pthread_barrier_t serverStart, serverEnd;
volatile bool serverEnabled = true;


void *serverThd(void *args)
{
#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();

    int rank = theTwoSided().getProcRank().toInt();

    CPU_SET(NUM_WORKERS, &cpuset);

    int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif

    // When enabled execute any active messages that arrive
    while (1) {
        pthread_barrier_wait(&serverStart);
        while (serverEnabled) {
            AmListenObjCodelet<NullMutex>* lcodelet;
            if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                lcodelet->execute();
            }
            Codelet* codelet;
            if ((codelet = serverDispatcher.progress()) != NULL) {
                codelet->execute();
                // Assume that codelet was an AmRequest and needs to be freed
                delete reinterpret_cast<AmRequest*>(codelet);
            }
        }
        pthread_barrier_wait(&serverEnd);
    }
    return NULL;
}

void *workerThd(void *args)
{
    const unsigned threadRank = *(unsigned*)args;

#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();

    int rank = theTwoSided().getProcRank().toInt();

    CPU_SET(threadRank, &cpuset);

    int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif

    while (1) {
        pthread_barrier_wait(&workersStart);
        utcs[threadRank]->process(threadRank);
        pthread_barrier_wait(&workersEnd);
    }

    return NULL;
}

void amBarrier()
{
    int epoch = pgrp->signalBarrier();
    while(!pgrp->testBarrier(epoch)) {
        AmListenObjCodelet<NullMutex>* codelet;
        if((codelet=theAm().amListeners[0]->progress()) != NULL) {
            codelet->execute();
        }
    }
}

void amBarrierThd()
{
    int epoch = pgrp->signalBarrier();
    while(!pgrp->testBarrier(epoch)) { }
}

static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}

typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_description;

static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_description *desc = (task_description*)_bigd;
    unsigned long seq_id[2];
    cell_t result;
    is_edge_param_t param;
    int is_edge_answer = 0;
    double t = 0;
    unsigned long i;
    int sscore;
    int maxLen;

    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
    {
        t = MPI_Wtime();
        affine_gap_align(
                sequences[seq_id[0]].c_str(), sequences[seq_id[0]].size(),
                sequences[seq_id[1]].c_str(), sequences[seq_id[1]].size(),
                &result, tbl[thd], del[thd], ins[thd]);
        param.AOL = 8;
        param.SIM = 4;
        param.OS = 3;
        is_edge_answer = is_edge(result,
                sequences[seq_id[0]].c_str(), sequences[seq_id[0]].size(),
                sequences[seq_id[1]].c_str(), sequences[seq_id[1]].size(),
                param, &sscore, &maxLen);
        ++stats[thd].align_counts;

        if (is_edge_answer) {
#if DEBUG
            cout << trank(thd)
                << ": aligned " << seq_id[0] << " " << seq_id[1]
                << ": (score,ndig,alen)=("
                << result.score << ","
                << result.ndig << ","
                << result.alen << ")"
                << ": edge? " << is_edge_answer << endl;
#endif
            ++stats[thd].edge_counts;
        }
        t = MPI_Wtime() - t;
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
    }
}

void createTaskFn(void *tsk, int tsk_size, TaskIndex tidx) {
    unsigned long seq_id[2];
    assert(tsk_size == sizeof(task_description));
    task_description *desc = (task_description*)tsk;
    k_combination(tidx, 2, seq_id);
    desc->id1 = seq_id[0];
    desc->id2 = seq_id[1];
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int provided;
    vector<string> all_argv;
    long file_size = -1;
    char *file_buffer = NULL;
    MPI_File fh;
    MPI_Status status;
    long seg_count = 0;
    size_t max_seq_len = 0;
    unsigned long nCk;

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

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, comm);
    tbl = new cell_t**[NUM_WORKERS];
    del = new int**[NUM_WORKERS];
    ins = new int**[NUM_WORKERS];
    utcs = new UniformTaskCollIter*[NUM_WORKERS];
    stats = new AlignStats[NUM_WORKERS];
    threadHandles = new pthread_t[NUM_WORKERS + NUM_SERVERS];
    threadRanks = new unsigned[NUM_WORKERS + NUM_SERVERS];
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        threadRanks[worker] = worker;
    }

    /* initialize dynamic code */
    init_map(SIGMA);

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
            int j;
            for (j=0; j<all_argv.size(); ++j) {
                printf("[%d] argv[%d]=%s\n", rank, j, all_argv[j].c_str());
            }
        }
        MPI_Barrier(comm);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() >= 4) {
                printf("too many arguments\n");
            }
            printf("usage: brute_force sequence_file [combinations_per_task]\n");
        }
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }

    /* process 0 open the file locally to determine its size */
    if (0 == rank) {
        FILE *file = fopen(all_argv[1].c_str(), "r");
        if (NULL == file) {
            printf("unable to open file on process 0\n");
            MPI_Abort(comm, 1);
        }
        if (0 != fseek(file, 0, SEEK_END)) {
            printf("unable to seek to end of file on process 0\n");
            MPI_Abort(comm, 1);
        }
        file_size = ftell(file);
        if (-1 == file_size) {
            printf("unable to get size of file on process 0\n");
            MPI_Abort(comm, 1);
        }
        if (0 != fclose(file)) {
            printf("unable to get close file on process 0\n");
            MPI_Abort(comm, 1);
        }
    }

    /* the file_size is broadcast to all */
    MPI_CHECK(MPI_Bcast(&file_size, 1, MPI_LONG, 0, comm));
    if (0 == trank(0)) {
        printf("file_size=%ld\n", file_size);
    }
    /* allocate a buffer for the file, of the entire size */
    /* TODO: this is not memory efficient since we allocate a buffer to read
     * the entire input file and then parse the buffer into a vector of
     * strings, essentially doubling the memory requirement */
    file_buffer = new char[file_size];

    /* all procs read the entire file */
    MPI_CHECK(MPI_File_open(comm, const_cast<char*>(all_argv[1].c_str()),
                MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN,
                MPI_INFO_NULL, &fh));
    MPI_CHECK(MPI_File_read_all(fh, file_buffer, file_size, MPI_CHAR, &status));
    MPI_CHECK(MPI_File_close(&fh));

    /* each process counts how many '>' characters are in the file_buffer */
    for (int i=0; i<file_size; ++i) {
        if (file_buffer[i] == '>') {
            ++seg_count;
        }
    }
#if DEBUG
    /* print the seg_count on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] seg_count=%ld\n", rank, seg_count);
        }
        MPI_Barrier(comm);
    }
#endif

    /* TODO declare these at the top */
    /* each process indexes the file_buffer */
    istrstream input_stream(const_cast<const char*>(file_buffer), file_size);
    string line;
    string sequence;
    sequences.reserve(seg_count);
    while (getline(input_stream, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                max_seq_len = max(max_seq_len, sequence.size());
            }
            sequence.clear();
            continue;
        }
        sequence += line;
    }
    /* add the last sequence in the file since we wouldn't encounter another
     * '>' character but rather an EOF */
    if (!sequence.empty()) {
        sequences.push_back(sequence);
        max_seq_len = max(max_seq_len, sequence.size());
    }
    sequence.clear();
#if DEBUG
    /* print the seg_count on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] sequences.size()=%ld sequences.capacity()=%ld\n",
                    rank, long(sequences.size()), long(sequences.capacity()));
        }
        MPI_Barrier(comm);
    }
    /* print the max_seq_len on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] max_seq_len=%ld\n", rank, long(max_seq_len));
        }
        MPI_Barrier(comm);
    }
#endif
#if DEBUG
    /* print the first sequence on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] sequences[0]=%s\n", rank, sequences[0].c_str());
        }
        MPI_Barrier(comm);
    }
    /* print the last sequence on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] sequences.back()=%s\n",
                    rank, sequences.back().c_str());
        }
        MPI_Barrier(comm);
    }
#endif
    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences.size(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences.size(), nCk);
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
    MPI_Barrier(comm);

    /* some more dynamic initialization */
    assert(NROW == 2);
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        tbl[worker] = alloc_tbl(NROW, max_seq_len);
        del[worker] = alloc_int(NROW, max_seq_len);
        ins[worker] = alloc_int(NROW, max_seq_len);
    }

    /* the tascel part */
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        int wrank = trank(worker);
        unsigned long lower_limit = wrank*tasks_per_worker;
        unsigned long upper_limit = lower_limit + tasks_per_worker;
        unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);
        /* if I'm the last worker, add the remainder of the tasks */
        if (wrank == nprocs*NUM_WORKERS-1) {
            upper_limit += remainder;
        }
        --upper_limit;
        UniformTaskCollIter*& utc = utcs[worker];
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description));
        utc = new UniformTaskCollIter(props, createTaskFn,
                lower_limit, upper_limit, worker);
    }

    amBarrier();

#if defined(THREADED)
#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();
    CPU_SET(0, &cpuset);
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif

    pthread_barrier_init(&workersStart, 0, NUM_WORKERS);
    pthread_barrier_init(&workersEnd, 0, NUM_WORKERS);
    pthread_barrier_init(&serverStart, 0, NUM_SERVERS + 1);
    pthread_barrier_init(&serverEnd, 0, NUM_SERVERS + 1);
    asm("mfence");
    for (unsigned i = 1; i < NUM_WORKERS; ++i) {
      pthread_create(&threadHandles[i], NULL, workerThd, &threadRanks[i]);
    }
    pthread_create(&threadHandles[NUM_WORKERS], NULL,
           serverThd, &threadRanks[NUM_WORKERS]);
#endif

#if defined(THREADED)
    serverEnabled = true;
    pthread_barrier_wait(&serverStart);
    MPI_Barrier(comm);
    pthread_barrier_wait(&workersStart);
#endif

    double mytimer = MPI_Wtime();
    utcs[0]->process(0);
    mytimer = MPI_Wtime() - mytimer;
    if (0 == trank(0)) {
        cout << "mytimer=" << mytimer << endl;
    }

#if defined(THREADED)
    pthread_barrier_wait(&workersEnd);
    amBarrierThd();

    serverEnabled = false;
    asm("mfence");
    pthread_barrier_wait(&serverEnd);

    amBarrier();
#endif

#if 0
    for (int i = 0; i < NUM_WORKERS; i++) {
        utcs[i]->restore();
    }
#endif

    amBarrier();
    MPI_Barrier(comm);

    AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
    MPI_Gather(stats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       rstats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank) {
        AlignStats totals;
        cout << " pid" << rstats[0].getHeader() << endl;      
        for(unsigned i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += rstats[i];
            cout << std::setw(4) << std::right << i << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT" << totals << endl;
    }
    delete [] rstats;
    rstats=NULL;

    StealingStats stt[NUM_WORKERS];
    for(unsigned i=0; i<NUM_WORKERS; i++) {
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
      for(unsigned i=0; i<nprocs*NUM_WORKERS; i++) {
	tstt += rstt[i];
	cout<<std::setw(4)<<std::right<<i<<rstt[i]<<endl;
      }
      cout<<"=============================================="<<endl;
      cout<<"TOT "<<tstt<<endl;
    }
    delete [] rstt;
    rstt=NULL;

    amBarrier();
    MPI_Barrier(comm);

    /* clean up */
    delete [] file_buffer;
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        free_tbl(tbl[worker], NROW);
        free_int(del[worker], NROW);
        free_int(ins[worker], NROW);
    }

    TascelConfig::finalize();
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
