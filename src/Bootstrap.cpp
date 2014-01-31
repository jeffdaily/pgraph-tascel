#include "config.h"

#include <cassert>

#include <mpi.h>

#include "Bootstrap.hpp"
#include "mpix.hpp"

#ifdef PAUSE_ON_ERROR
#include "SigSegvHandler.hpp"
#endif


namespace pgraph {

MPI_Comm comm = MPI_COMM_NULL;
int rank = 0;
int nprocs = 0;

void initialize(int *argc, char ***argv)
{
#ifdef PAUSE_ON_ERROR
    TrapSigSegv();
#endif

    /* initialize MPI */
#if defined(THREADED)
    {
        int provided;
        MPI_CHECK(MPI_Init_thread(argc, argv,
                    MPI_THREAD_MULTIPLE, &provided));
        assert(provided == MPI_THREAD_MULTIPLE);
    }
#else
    MPI_CHECK(MPI_Init(argc, argv));
#endif

    /* get rank and nprocs */
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));
    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &nprocs));
}


void finalize()
{
    MPI_CHECK(MPI_Comm_free(&comm));
    MPI_CHECK(MPI_Finalize());
}

};

