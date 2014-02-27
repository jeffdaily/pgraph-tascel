#include "config.h"

#include <cassert>

#include <mpi.h>

#include "Bootstrap.hpp"
#include "mpix.hpp"

#define PAUSE_ON_ERROR 1
#ifdef PAUSE_ON_ERROR
#include "SigSegvHandler.hpp"
#endif

#if HAVE_ARMCI
/* ANL's armci does not extern "C" inside the header */
extern "C" {
#include <armci.h>
}
#endif

namespace pgraph {

MPI_Comm comm = MPI_COMM_NULL;
int rank = 0;
int nprocs = 0;
int *argc_ = NULL;
char ***argv_ = NULL;
bool is_armci_initialized = false;

void initialize(int *argc, char ***argv)
{
#ifdef PAUSE_ON_ERROR
    TrapSigSegv();
#endif

    argc_ = argc;
    argv_ = argv;
    is_armci_initialized = false;

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
    finalize_armci();
    MPI_CHECK(MPI_Comm_free(&comm));
    MPI_CHECK(MPI_Finalize());
}


void initialize_armci()
{
#if HAVE_ARMCI
    if (!is_armci_initialized) {
        ARMCI_Init_args(argc_,argv_);
        is_armci_initialized = true;
    }
#endif
}


void finalize_armci()
{
#if HAVE_ARMCI
    if (is_armci_initialized) {
        ARMCI_Finalize();
    }
#endif
}

};

