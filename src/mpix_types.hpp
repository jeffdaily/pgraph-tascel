/**
 * @file mpix_types.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * MPI_Datatype creation functions for pgraph structures.
 */
#ifndef _MPIX_TYPES_HPP_
#define _MPIX_TYPES_HPP_

#include <mpi.h>

#include <tascel/Counter.h>
#include <tascel/StealingStats.h>
#include <tascel/Timer.h>

#include "AlignStats.hpp"
#include "DbStats.hpp"
#include "DupStats.hpp"
#include "Stats.hpp"
#include "Suffix.hpp"
#include "SuffixArrayStats.hpp"
#include "TreeStats.hpp"

#include "mpix.hpp"
#include "mpix_helper.hpp"

using ::pgraph::AlignStats;
using ::pgraph::DbStats;
using ::pgraph::DupStats;
using ::pgraph::Stats;
using ::pgraph::Suffix;
using ::pgraph::SuffixArrayStats;
using ::pgraph::TreeStats;

using ::tascel::Counter;
using ::tascel::StealingStats;
using ::tascel::Timer;

namespace mpix {

extern MPI_Datatype mpi_datatype_Counter;
extern MPI_Datatype mpi_datatype_Timer;
extern MPI_Datatype mpi_datatype_StealingStats;
extern MPI_Datatype mpi_datatype_Stats;
extern MPI_Datatype mpi_datatype_SuffixArrayStats;
extern MPI_Datatype mpi_datatype_TreeStats;
extern MPI_Datatype mpi_datatype_DbStats;
extern MPI_Datatype mpi_datatype_DupStats;
extern MPI_Datatype mpi_datatype_AlignStats;
extern MPI_Datatype mpi_datatype_Suffix;

MPI_Datatype get_mpi_datatype(Counter object);
MPI_Datatype get_mpi_datatype(Timer object);
MPI_Datatype get_mpi_datatype(StealingStats object);
MPI_Datatype get_mpi_datatype(Stats object);
MPI_Datatype get_mpi_datatype(SuffixArrayStats object);
MPI_Datatype get_mpi_datatype(TreeStats object);
MPI_Datatype get_mpi_datatype(DbStats object);
MPI_Datatype get_mpi_datatype(DupStats object);
MPI_Datatype get_mpi_datatype(AlignStats object);
MPI_Datatype get_mpi_datatype(Suffix object);

void init_types();
void free_types();

MPIX_DECL_ALL(Stats)
MPIX_DECL_ALL(SuffixArrayStats)
MPIX_DECL_ALL(TreeStats)
MPIX_DECL_ALL(DbStats)
MPIX_DECL_ALL(DupStats)
MPIX_DECL_ALL(AlignStats)
MPIX_DECL_ALL(StealingStats)

}

#endif /* _MPIX_TYPES_HPP_ */
