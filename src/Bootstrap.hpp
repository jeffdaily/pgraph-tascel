#ifndef _PGRAPH_BOOTSTRAP_H_
#define _PGRAPH_BOOTSTRAP_H_

#include <mpi.h>

namespace pgraph {
    extern MPI_Comm comm;
    extern int rank;
    extern int nprocs;
    void initialize(int *argc, char ***argv);
    void finalize();
};

#endif /* _PGRAPH_BOOTSTRAP_H_ */
