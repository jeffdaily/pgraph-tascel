#ifndef _NXTVAL_H_
#define _NXTVAL_H_

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include <pthread.h>

#define TAG_REQUEST 698825
#define TAG_RESPONSE 698826

void* nxtval_server(void *ignore)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    long long nxtval = *((long long*)ignore);
    long long recval = 0;
    MPI_Status status;
    free(ignore);
    while (1) {
        MPI_Recv(&recval, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, TAG_REQUEST, comm, &status);
        if (recval <= 0) {
            printf("nxtval_server got break\n");
            break;
        }
        MPI_Send(&nxtval, 1, MPI_LONG_LONG, status.MPI_SOURCE, TAG_RESPONSE, comm);
        nxtval += recval;
    }
    return NULL;
}


typedef struct NXTVAL {
    MPI_Comm comm;
    MPI_Comm comm2;
    int rank;
    int size;
    int provided;
    pthread_t thread;
    long long start;
    long long stop;
} NXTVAL_t;


static NXTVAL_t NXTVAL_init(long long start, long long stop)
{
    NXTVAL_t nxt;

    nxt.comm = MPI_COMM_WORLD;
    MPI_Comm_rank(nxt.comm, &nxt.rank);
    MPI_Comm_size(nxt.comm, &nxt.size);
    nxt.start = start;
    nxt.stop = stop;

    MPI_Barrier(nxt.comm);

    if (0 == nxt.rank) {
        printf("NXTVAL_init (size=%d)\n", nxt.size);
        fflush(stdout);
    }

    MPI_Query_thread(&nxt.provided);

    if (nxt.size > 1) {
        if (MPI_THREAD_MULTIPLE == nxt.provided) {
            if (0 == nxt.rank) {
                int rc;
                long long *arg = (long long*)malloc(sizeof(long long));
                *arg = start;
                printf("starting nxtval_server as thread\n");
                fflush(stdout);
                rc = pthread_create(&nxt.thread, NULL, nxtval_server, arg);
                if (rc) {
                    printf("ERROR return code from pthread_create() is %d\n", rc);
                    fflush(stdout);
                    MPI_Abort(nxt.comm, -1);
                    exit(EXIT_FAILURE);
                }
            }
        }
        else {
            if (0 == nxt.rank) {
                long long *arg = (long long*)malloc(sizeof(long long));
                *arg = start;
                printf("starting nxtval_server as rank 0\n");
                fflush(stdout);
                MPI_Comm_split(nxt.comm, 0, nxt.rank, &nxt.comm2);
                (void)nxtval_server(arg);
            }
            else {
                MPI_Comm_split(nxt.comm, 1, nxt.rank, &nxt.comm2);
            }
        }
    }
    else {
        printf("NXTVAL_init not started due to single rank\n");
        fflush(stdout);
    }

    return nxt;
}


static void NXTVAL_stop(NXTVAL_t nxt)
{
    printf("%d: stopping nxtval\n", nxt.rank);
    fflush(stdout);

    if (nxt.size > 1) {
        if (MPI_THREAD_MULTIPLE == nxt.provided) {
            MPI_Barrier(nxt.comm);
            if (0 == nxt.rank) {
                /* rank 0 stops the thread */
                void *res;
                int rc;
                long long zero = 0;

                printf("rank 0 stopping thread\n");
                fflush(stdout);

                MPI_Send(&zero, 1, MPI_LONG_LONG, 0, TAG_REQUEST, nxt.comm);
                rc = pthread_join(nxt.thread, &res);
                if (rc) {
                    printf("ERROR return code from pthread_join() is %d\n", rc);
                    fflush(stdout);
                    MPI_Abort(nxt.comm, -1);
                    exit(EXIT_FAILURE);
                }
            }
        }
        else {
            MPI_Barrier(nxt.comm2);
            if (1 == nxt.rank) {
                /* rank 1 stops rank 0 */
                long long zero = 0;

                printf("rank 1 stopping rank 0\n");
                fflush(stdout);

                MPI_Send(&zero, 1, MPI_LONG_LONG, 0, TAG_REQUEST, nxt.comm);
            }
            MPI_Comm_free(&nxt.comm2);
        }

        MPI_Barrier(nxt.comm);
    }
    else {
        printf("NXTVAL_stop not needed, single rank\n");
        fflush(stdout);
    }
}


static long long NXTVAL_get(NXTVAL_t nxt)
{
    if (nxt.size > 1) {
        if (MPI_THREAD_MULTIPLE != nxt.provided && 0 == nxt.rank) {
            return nxt.stop;
        }
        else {
            MPI_Status status;
            long long one=1;
            long long oldval;
            MPI_Send(&one, 1, MPI_LONG_LONG, 0, TAG_REQUEST, nxt.comm);
            MPI_Recv(&oldval, 1, MPI_LONG_LONG, 0, TAG_RESPONSE, nxt.comm, &status);
            return oldval;
        }
    }
    else {
        static long long nxtval = nxt.start;
        long long oldval = nxtval;
        nxtval += 1;
        return oldval;
    }
}

#endif /* _NXTVAL_H_ */

