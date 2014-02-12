/**
 * @author jeff.daily@pnnl.gov
 *
 * Alignment of an input dataset using work stealing. Each MPI task reads the
 * input file.
 */
#include "config.h"

#include <mpi.h>
extern "C" {
#include <armci.h>
}

#include <fcntl.h>
#include <sys/errno.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>
#include <semaphore.h>

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

#include "alignment.hpp"
#include "AlignStats.hpp"
#include "Bootstrap.hpp"
#include "combinations.h"
#include "constants.h"
#include "EdgeResult.hpp"
#include "mpix.hpp"
#include "Parameters.hpp"

#define USE_SEQ_DB 1
#if USE_SEQ_DB
#include "SequenceDatabaseArmci.hpp"
#endif

using namespace std;
using namespace pgraph;

#define SEP ","
#define ALL_RESULTS 1
#define GB (1073741824)

typedef struct {
    const char *str;
    size_t size;
} sequence_t;

//static string _get_edges_filename(int rank);
#if !USE_SEQ_DB
static unsigned long _mytime();
static string _generate_shm_name();
static void* _shm_create(MPI_Comm comm, size_t size, string &name);
static void _pack_and_index_fasta(char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t &max_seq_size);
static void _index_fasta(const char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t &max_seq_size);
#endif
static size_t parse_memory_budget(const string& value);


int main(int argc, char **argv)
{
    MPI_Comm comm_world = MPI_COMM_NULL;
    int rank_world = 0;
    int size_world = 0;
    AlignStats stats;
    vector<EdgeResult> edge_results;
    Parameters parameters;
    vector<string> all_argv;
    unsigned long nCk = 0;
    size_t sequence_count = 0;
    int retval = 0;
    size_t max_seq_size = 0;
#if !USE_SEQ_DB
    MPI_Comm comm_node = MPI_COMM_NULL;
    MPI_Comm comm_master = MPI_COMM_NULL;
    int rank_node = 0;
    int rank_master = 0;
    int size_node = 0;
    int size_master = 0;
    MPI_Offset file_size = -1;
    char *file_buffer = NULL;
    vector<sequence_t> sequences;
    string shm_name;
#else
    SequenceDatabaseArmci *sequences = NULL;
#endif
#if USE_SSW
#else
    cell_t **tbl = NULL;
    int **ins = NULL;
    int **del = NULL;
#endif

    /* initialize MPI */
    MPI_CHECK_C(MPI_Init(&argc, &argv));
    MPI_CHECK_C(MPI_Comm_dup(MPI_COMM_WORLD, &comm_world));
    MPI_CHECK_C(MPI_Comm_rank(comm_world, &rank_world));
    MPI_CHECK_C(MPI_Comm_size(comm_world, &size_world));
    
    /* initialize ARMCI */
    ARMCI_Init_args(&argc, &argv);

#if !USE_SEQ_DB
    /* split comm_world based on hostid */
    MPI_CHECK_C(MPI_Comm_split(comm_world, gethostid(), rank_world, &comm_node));
    MPI_CHECK_C(MPI_Comm_rank(comm_node, &rank_node));
    MPI_CHECK_C(MPI_Comm_size(comm_node, &size_node));

    /* split comm_world based on node master */
    if (0 == rank_node) {
        MPI_CHECK_C(MPI_Comm_split(comm_world, 0, rank_world, &comm_master));
        MPI_CHECK_C(MPI_Comm_rank(comm_master, &rank_master));
        MPI_CHECK_C(MPI_Comm_size(comm_master, &size_master));
    }
    else {
        MPI_CHECK_C(MPI_Comm_split(comm_world, 1, rank_world, &comm_master));
        MPI_CHECK_C(MPI_Comm_free(&comm_master));
        comm_master = MPI_COMM_NULL;
        rank_master = -1;
        size_master = -1;
    }

#if DEBUG
    mpix_print_sync("rank_world", rank_world, comm_world);
    mpix_print_sync("rank_node",  rank_node, comm_world);
    mpix_print_sync("rank_master",rank_master, comm_world);
    mpix_print_sync("size_world", size_world, comm_world);
    mpix_print_sync("size_node",  size_node, comm_world);
    mpix_print_sync("size_master",size_master, comm_world);
#endif
#endif

    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank_world) {
        MPI_CHECK_C(MPI_Bcast(&argc, 1, MPI_INT, 0, comm_world));
        for (int i=0; i<argc; ++i) {
            int length = strlen(argv[i])+1;
            MPI_CHECK_C(MPI_Bcast(&length, 1, MPI_INT, 0, comm_world));
            MPI_CHECK_C(MPI_Bcast(argv[i], length, MPI_CHAR, 0, comm_world));
            all_argv.push_back(argv[i]);
        }
    } else {
        int all_argc;
        MPI_CHECK_C(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm_world));
        for (int i=0; i<all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK_C(MPI_Bcast(&length, 1, MPI_INT, 0, comm_world));
            MPI_CHECK_C(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm_world));
            all_argv.push_back(buffer);
        }
    }

#if DEBUG
    /* print the command line arguments */
    for (int i=0; i<size_world; ++i) {
        if (i == rank_world) {
            for (size_t j=0; j<all_argv.size(); ++j) {
                printf("[%d] argv[%zd]=%s\n", rank_world, j, all_argv[j].c_str());
            }
        }
        MPI_Barrier(comm_world);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 2 || all_argv.size() >= 5) {
        if (0 == rank_world) {
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
        parameters.parse(all_argv[3].c_str(), comm_world);
    } 
    if (0 == rank_world) {
        printf("----------------------------------------------\n");
        printf("%-20s: %d\n", "slide size", parameters.window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters.exact_match_length);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters.AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters.SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters.OS);
        printf("----------------------------------------------\n");
    }

#if USE_SEQ_DB
    sequences = new SequenceDatabaseArmci(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm_world, 1, '\0');
    max_seq_size = sequences->get_max_length();
#else
    /* rank_world 0 on each node will open the file into shared memory */
    if (0 == rank_node) {
        file_size = mpix_get_file_size(all_argv[1], comm_master);
    }
    mpix_bcast(file_size, 0, comm_node);
    file_buffer = (char*)_shm_create(comm_node, file_size, shm_name);

    /* rank_world 0 on each node will open the file into shared memory */
    if (0 == rank_node) {
        mpix_read_file(all_argv[1], file_buffer, file_size, GB, comm_master);
        _pack_and_index_fasta(file_buffer, file_size, '\0', sequences, max_seq_size);
        retval = msync(file_buffer, file_size, MS_SYNC|MS_INVALIDATE);
        if (-1 == retval) {
            perror("msync");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else {
        while (file_buffer[0] != '>') {
            //cout << "waiting for shm file" << endl;
            sleep(1);
        }
        _index_fasta(file_buffer, file_size, '\0', sequences, max_seq_size);
    }
    MPI_Barrier(comm_world);
#endif

#if DEBUG
    mpix_print_sync("max_seq_size", max_seq_size, comm_world);
#endif

#if USE_SSW
#else
    {
        size_t max_seq_size_p1 = max_seq_size + 1;
        tbl = new cell_t*[2];
        tbl[0] = new cell_t[max_seq_size_p1];
        tbl[1] = new cell_t[max_seq_size_p1];
        del = new int*[2];
        del[0] = new int[max_seq_size_p1];
        del[1] = new int[max_seq_size_p1];
        ins = new int*[2];
        ins[0] = new int[max_seq_size_p1];
        ins[1] = new int[max_seq_size_p1];
#if USE_MEMSET
        (void)memset(tbl[0], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(tbl[1], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(del[0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(del[1], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[1], 0, sizeof(int)*max_seq_size_p1);
#endif
    }
#endif

    /* how many combinations of sequences are there? */
#if USE_SEQ_DB
    sequence_count = sequences->get_global_count();
#else
    sequence_count = sequences.size();
#endif
    nCk = binomial_coefficient(sequence_count, 2);
    if (0 == rank_world) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequence_count, nCk);
    }

    double selectivity = 1.0;
    if (all_argv.size() == 3) {
      selectivity  = fabs(min(1.0,atof(all_argv[2].c_str())));
    }
    long nalignments = (long)(0.5+selectivity*nCk);
    if (0 == rank_world) {
        printf("selectivity=%lf\n", selectivity);
        printf("nalignments=%ld\n", nalignments);
    }
    MPI_Barrier(comm_world);

    /* allocate task counter */
    long task_id = rank_world;
    long bytes = 0;
    long **counter = NULL;

    counter = new long*[size_world];
    bytes = (0 == rank_world) ? sizeof(long) : 0;
    retval = ARMCI_Malloc((void**)counter, bytes);
    assert(0 == retval);
    ARMCI_Barrier();
    if (0 == rank_world) {
        counter[0][0] = size_world; /* init counter */
    }
    ARMCI_Barrier();

    double mytimer = MPI_Wtime();
    double mytime_combo = 0.0;
    double mytime_fetch = 0.0;
    double mytime_misc = 0.0;
    double mytime_total = 0.0;
    // TODO TAKS COUNTER AND ALIGNMENT LOOP
    while (task_id < nalignments) {
        unsigned long seq_id[2];
        double timer_total = MPI_Wtime();
        double timer;
        timer = MPI_Wtime();
        k_combination2(task_id, seq_id);
        timer = MPI_Wtime() - timer;
        mytime_combo += timer;
        stats.kcomb_times_tot += timer;
        timer = MPI_Wtime();
        assert(seq_id[0] < sequence_count);
        assert(seq_id[1] < sequence_count);
#if !USE_SEQ_DB
        assert(sequences[seq_id[0]].str);
        assert(sequences[seq_id[1]].str);
        assert(sequences[seq_id[0]].size);
        assert(sequences[seq_id[1]].size);
        assert(sequences[seq_id[0]].size <= max_seq_size);
        assert(sequences[seq_id[1]].size <= max_seq_size);
#endif
        unsigned long s1Len = (*sequences)[seq_id[0]].get_sequence_length();
        unsigned long s2Len = (*sequences)[seq_id[1]].get_sequence_length();
#ifdef LENGTH_FILTER
        int cutOff = parameters.AOL * parameters.SIM;
        if ((s1Len <= s2Len && (100 * s1Len < cutOff * s2Len))
                || (s2Len < s1Len && (100 * s2Len < cutOff * s1Len))) {
            stats.work_skipped += s1Len * s2Len;
            ++stats.align_skipped;
        }
        else
#endif
        {
            cell_t result;
            double t = MPI_Wtime();
            int sscore;
            size_t maxLen;
            stats.work += s1Len * s2Len;
#if USE_SSW
#if USE_SEQ_DB
            sequences->align_ssw(seq_id[0], seq_id[1], result.score, result.matches, result.length, parameters.open, parameters.gap);

#else
            result = pgraph::affine_gap_align_blosum_ssw(
                    sequences[seq_id[0]].str, sequences[seq_id[0]].size,
                    sequences[seq_id[1]].str, sequences[seq_id[1]].size,
                    parameters.open, parameters.gap);
#endif
#else
#if USE_SEQ_DB
            sequences->align(seq_id[0], seq_id[1], result.score, result.matches, result.length, parameters.open, parameters.gap);
#else
            result = pgraph::affine_gap_align_blosum(
                    sequences[seq_id[0]].str, sequences[seq_id[0]].size,
                    sequences[seq_id[1]].str, sequences[seq_id[1]].size,
                    tbl, del, ins, parameters.open, parameters.gap);
#endif
#endif
#if USE_SEQ_DB
            bool is_edge_answer = sequences->is_edge(
                    seq_id[0],
                    seq_id[1],
                    result.score, result.matches, result.length,
                    parameters.AOL, parameters.SIM, parameters.OS,
                    sscore, maxLen);
#else
            bool is_edge_answer = pgraph::is_edge_blosum(result,
                    sequences[seq_id[0]].str, sequences[seq_id[0]].size,
                    sequences[seq_id[1]].str, sequences[seq_id[1]].size,
                    parameters.AOL, parameters.SIM, parameters.OS,
                    sscore, maxLen);
#endif

            ++stats.align_counts;

            if (is_edge_answer || ALL_RESULTS)
            {
#if DEBUG
                cout << rank_world
                    << ": aligned " << seq_id[0] << " " << seq_id[1]
                    << ": (score,matches,length)=("
                    << result.score << ","
                    << result.matches << ","
                    << result.length << ")"
                    << ": edge? " << is_edge_answer << endl;
#endif
                edge_results.push_back(EdgeResult(
                            seq_id[0], seq_id[1],
#if 0
                            1.0*result.length/max_len,
                            1.0*result.matches/result.length,
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
                    ++stats.edge_counts;
                }
            }
            t = MPI_Wtime() - t;
            stats.align_times_tot += t;
            stats.calc_min(t);
            stats.calc_max(t);
        }

        // next task
        timer = MPI_Wtime()- timer;
        mytime_misc += timer;
        timer = MPI_Wtime();
        retval = ARMCI_Rmw(ARMCI_FETCH_AND_ADD_LONG, &task_id, counter[0], 1, 0);
        timer = MPI_Wtime() - timer;
        mytime_fetch += timer;
        assert(0 == retval);
        timer_total = MPI_Wtime() - timer_total;
        mytime_total += timer_total;
        stats.total_times_tot += timer_total;
    }

    //if (0 == rank_world) {
        retval = ARMCI_Free(counter[rank_world]);
        delete [] counter;
    //}

    mytimer = MPI_Wtime() - mytimer;

    if (0 == rank_world) {
        cout << "mytimer=" << mytimer << endl;
        cout << "mytime_combo=" << mytime_combo << endl;
        cout << "mytime_fetch=" << mytime_fetch << endl;
        cout << "mytime_misc=" << mytime_misc << endl;
        cout << "mytime_total=" << mytime_total << endl;
    }

    MPI_Barrier(comm_world);

    AlignStats * rstats = new AlignStats[size_world];
    MPI_Gather(&stats, sizeof(AlignStats), MPI_CHAR, 
	       rstats, sizeof(AlignStats), MPI_CHAR, 
	       0, comm_world);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank_world) {
        AlignStats totals;
        cout << " pid" << rstats[0].getHeader() << endl;      
        for(int i=0; i<size_world; i++) {
            totals += rstats[i];
            cout << std::setw(4) << std::right << i << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT" << totals << endl;
    }
    delete [] rstats;
    rstats=NULL;

    double * combo_stats = new double[size_world];
    MPI_Gather(&mytime_combo, 1, MPI_DOUBLE,
            combo_stats, 1, MPI_DOUBLE,
            0, comm_world);
    if (0 == rank_world) {
        double totals = 0.0;
        cout << " pid k_combination_time" << endl;
        for (int i=0; i<size_world; ++i) {
            totals += combo_stats[i];
            //cout << std::setw(4) << std::right << i << " " << combo_stats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT " << totals << " AVG " << totals/size_world << endl;
    }
    delete [] combo_stats;

    double * fetch_stats = new double[size_world];
    MPI_Gather(&mytime_fetch, 1, MPI_DOUBLE,
            fetch_stats, 1, MPI_DOUBLE,
            0, comm_world);
    if (0 == rank_world) {
        double totals = 0.0;
        cout << " pid fetch_time" << endl;
        for (int i=0; i<size_world; ++i) {
            totals += fetch_stats[i];
            //cout << std::setw(4) << std::right << i << " " << fetch_stats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT " << totals << " AVG " << totals/size_world << endl;
    }
    delete [] fetch_stats;

    MPI_Barrier(comm_world);

#if !USE_SEQ_DB
    /* unmap the memory */
    retval = munmap(file_buffer, file_size);
    if (-1 == retval) {
        perror("munmap");
        //MPI_Abort(MPI_COMM_WORLD, retval);
    }

    /* remove the shared memory object */
    if (0 == rank_world) {
        retval = shm_unlink(shm_name.c_str());
        if (-1 == retval) {
            perror("shm_unlink");
            //MPI_Abort(MPI_COMM_WORLD, retval);
        }
    }
#endif

#if USE_SSW
#else
    /* clean up DP tables */
    delete [] tbl[0];
    delete [] tbl[1];
    delete [] del[0];
    delete [] del[1];
    delete [] ins[0];
    delete [] ins[1];
    delete [] tbl;
    delete [] del;
    delete [] ins;
#endif
#if USE_SEQ_DB
    delete sequences;
#endif

    ARMCI_Finalize();

#if !USE_SEQ_DB
    if (MPI_COMM_NULL != comm_master) {
        MPI_Comm_free(&comm_master);
    }
    MPI_Comm_free(&comm_node);
#endif
    pgraph::finalize();

    return 0;
}


#if 0
static string _get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}
#endif


#if !USE_SEQ_DB
static void* _shm_create(MPI_Comm comm, size_t size, string &name)
{
    void *mapped = NULL;
    int fd = 0;
    int retval = 0;
    int rank = 0;

    MPI_CHECK_C(MPI_Comm_rank(comm, &rank));

    if (0 == rank) {
        name = _generate_shm_name();

        /* create shared memory segment */
        fd = shm_open(name.c_str(), O_CREAT|O_EXCL|O_RDWR, S_IRUSR|S_IWUSR);
        if (-1 == fd) {
            perror("_shm_create: shm_open");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* set the size of my shared memory object */
        retval = ftruncate(fd, size);
        if (-1 == retval) {
            perror("_shm_create: ftruncate");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* map into local address space */
        mapped = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        if (MAP_FAILED == mapped) {
            perror("_shm_map: mmap");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* close file descriptor */
        retval = close(fd);
        if (-1 == retval) {
            perror("_shm_create: close");
            MPI_Abort(MPI_COMM_WORLD, retval);
        }

        /* broadcast the name */
        mpix_bcast(name, 0, comm);
    }
    else {
        /* broadcast the name */
        mpix_bcast(name, 0, comm);

        /* attach to shared memory segment */
        fd = shm_open(name.c_str(), O_RDWR, S_IRUSR|S_IWUSR);
        if (-1 == fd) {
            perror("_shm_create (attach): shm_open");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* map into local address space */
        mapped = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        if (MAP_FAILED == mapped) {
            perror("_shm_map: mmap");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* close file descriptor */
        retval = close(fd);
        if (-1 == retval) {
            perror("_shm_create (attach): close");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    return mapped;
}


static unsigned long _mytime()
{
    struct timeval tv;
    int retval = 0;

    retval = gettimeofday(&tv, NULL);
    assert(0 == retval);

    return ((unsigned long)tv.tv_sec)*1000000UL + ((unsigned long)tv.tv_usec);
}


static string _generate_shm_name()
{
    ostringstream os;
    os << "/pgraph" << _mytime() << getpid();
    return os.str();
}


static void _pack_and_index_fasta(char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t &max_seq_size)
{
    size_t r = 0;           /* read index */
    size_t w = 0;           /* write index */
    size_t last_gt = 0;     /* index of last '>' seen */
    size_t last_hash = 0;   /* index of last '#' inserted */
    sequence_t sequence;

    /* We skip all newlines unless it is the newline that terminates the
     * sequence ID. We replace the ID-terminating newline with '#'. We replace
     * the sequence terminating newline with the given delimiter. */

    assert(buffer[0] == '>');

    while (r<size) {
        if (buffer[r] == '>') {
            last_gt = w;
            while (r<size && buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            last_hash = w;
            buffer[w++] = '#';
            r++;
        }
        else {
            while (r<size && buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            /* peek at next character, either EOF or '>' */
            if (r<size) {
                if (r+1 >= size || buffer[r+1] == '>') {
                    buffer[w++] = delimiter;
                    size_t sequence_offset = last_hash + 1;
                    size_t sequence_length = w - last_hash - 1;
                    if (delimiter == '\0') {
                        sequence_length -= 1;
                    }
                    sequence.str = &buffer[sequence_offset];
                    sequence.size = sequence_length;
                    sequences.push_back(sequence);
                    max_seq_size = sequence.size > max_seq_size ?
                        sequence.size : max_seq_size;
                }
            }
            /* or file wasn't terminated with a newline */
            else {
                assert(0);
            }
            r++;
        }
    }
}


static void _index_fasta(const char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t &max_seq_size)
{
    size_t r = 0;           /* read index */
    size_t last_gt = 0;     /* index of last '>' seen */
    size_t last_hash = 0;   /* index of last '#' inserted */
    sequence_t sequence;

    /* We skip all newlines unless it is the newline that terminates the
     * sequence ID. We replace the ID-terminating newline with '#'. We replace
     * the sequence terminating newline with the given delimiter. */

    assert(buffer[0] == '>');

    while (r<size) {
        if (buffer[r] == '>') {
            if (r > 0) {
                size_t sequence_offset = last_hash + 1;
                size_t sequence_length = r - last_hash - 1;
                if (delimiter == '\0') {
                    sequence_length -= 1;
                }
                sequence.str = &buffer[sequence_offset];
                sequence.size = sequence_length;
                sequences.push_back(sequence);
                max_seq_size = sequence.size > max_seq_size ?
                    sequence.size : max_seq_size;
            }
            last_gt = r;
        }
        else if (buffer[r] == '#') {
            last_hash = r;
        }
        ++r;
    }
    {
        size_t sequence_offset = last_hash + 1;
        size_t sequence_length = r - last_hash - 1;
        if (delimiter == '\0') {
            sequence_length -= 1;
        }
        sequence.str = &buffer[sequence_offset];
        sequence.size = sequence_length;
        sequences.push_back(sequence);
        max_seq_size = sequence.size > max_seq_size ?
            sequence.size : max_seq_size;
    }
}
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

