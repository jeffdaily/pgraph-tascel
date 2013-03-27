ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS =
# this flag is supported by g++ and icpc
AM_CXXFLAGS = -Wno-deprecated
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD =

AM_CPPFLAGS += $(TASCEL_CPPFLAGS)
#AM_CPPFLAGS += $(GA_CPPFLAGS)
AM_CPPFLAGS += $(ARMCI_CPPFLAGS)
AM_CPPFLAGS += $(PT_MPI_CPPFLAGS)
AM_CPPFLAGS += $(GMP_CPPFLAGS)
AM_CPPFLAGS += -I$(top_srcdir)/src
AM_CPPFLAGS += -I$(top_srcdir)/src/suffixtree

AM_LDFLAGS += $(TASCEL_LDFLAGS)
#AM_LDFLAGS += $(GA_LDFLAGS)
AM_LDFLAGS += $(ARMCI_LDFLAGS)
AM_LDFLAGS += $(PT_MPI_LDFLAGS)
AM_LDFLAGS += $(GMP_LDFLAGS)

LDADD += libpgraph.la
LDADD += $(TASCEL_LIBS)
#LDADD += $(GA_LIBS)
LDADD += $(ARMCI_LIBS)
LDADD += $(PT_MPI_LIBS)
LDADD += $(GMP_LIBS)
#LDADD += $(GA_FLIBS)

noinst_LTLIBRARIES =
noinst_LTLIBRARIES += libpgraph.la

libpgraph_la_SOURCES =
libpgraph_la_SOURCES += src/AlignStats.h
libpgraph_la_SOURCES += src/blosum/blosum40.h
libpgraph_la_SOURCES += src/blosum/blosum45.h
libpgraph_la_SOURCES += src/blosum/blosum62.h
libpgraph_la_SOURCES += src/blosum/blosum75.h
libpgraph_la_SOURCES += src/blosum/blosum80.h
libpgraph_la_SOURCES += src/blosum/blosum90.h
libpgraph_la_SOURCES += src/combinations.c
libpgraph_la_SOURCES += src/combinations.h
libpgraph_la_SOURCES += src/dynamic.c
libpgraph_la_SOURCES += src/dynamic.h
libpgraph_la_SOURCES += src/io.cc
libpgraph_la_SOURCES += src/io.h
libpgraph_la_SOURCES += src/mpix.cc
libpgraph_la_SOURCES += src/mpix.h
libpgraph_la_SOURCES += src/pthread_fixes.h
libpgraph_la_SOURCES += src/Sequence.cc
libpgraph_la_SOURCES += src/Sequence.h
libpgraph_la_SOURCES += src/SequenceDatabase.cc
libpgraph_la_SOURCES += src/SequenceDatabase.h
libpgraph_la_SOURCES += src/suffixtree/bucket.c
libpgraph_la_SOURCES += src/suffixtree/bucket.h
libpgraph_la_SOURCES += src/suffixtree/cfg.c
libpgraph_la_SOURCES += src/suffixtree/cfg.h
libpgraph_la_SOURCES += src/suffixtree/dynamic.c
libpgraph_la_SOURCES += src/suffixtree/dynamic.h
libpgraph_la_SOURCES += src/suffixtree/elib.c
libpgraph_la_SOURCES += src/suffixtree/elib.h
libpgraph_la_SOURCES += src/suffixtree/lib.c
libpgraph_la_SOURCES += src/suffixtree/lib.h
libpgraph_la_SOURCES += src/suffixtree/loadseq.c
libpgraph_la_SOURCES += src/suffixtree/loadseq.h
libpgraph_la_SOURCES += src/suffixtree/pairs.c
libpgraph_la_SOURCES += src/suffixtree/pairs.h
libpgraph_la_SOURCES += src/suffixtree/search.c
libpgraph_la_SOURCES += src/suffixtree/search.h
libpgraph_la_SOURCES += src/suffixtree/stree.c
libpgraph_la_SOURCES += src/suffixtree/stree.h
libpgraph_la_SOURCES += src/suffixtree/type.h
libpgraph_la_SOURCES += src/suffixtree/uFind.c
libpgraph_la_SOURCES += src/suffixtree/uFind.h
libpgraph_la_SOURCES += src/suffixtree/wlib.c
libpgraph_la_SOURCES += src/suffixtree/wlib.h
libpgraph_la_SOURCES += src/timer.h

bin_PROGRAMS =
bin_PROGRAMS += apps/align
bin_PROGRAMS += apps/ccd

noinst_PROGRAMS =
noinst_PROGRAMS += sandbox/blast
noinst_PROGRAMS += sandbox/brute_force
if HAVE_ARMCI
noinst_PROGRAMS += sandbox/brute_force_armci
endif
noinst_PROGRAMS += sandbox/brute_force_binning
noinst_PROGRAMS += sandbox/brute_force_costs
noinst_PROGRAMS += sandbox/brute_force_multiple
noinst_PROGRAMS += sandbox/brute_force_pairs
noinst_PROGRAMS += sandbox/brute_force_rr
noinst_PROGRAMS += sandbox/brute_force_rr_sorted
noinst_PROGRAMS += sandbox/brute_force_rr_sorted_seq
noinst_PROGRAMS += sandbox/brute_force_serial
noinst_PROGRAMS += sandbox/brute_force_sorted
noinst_PROGRAMS += sandbox/brute_force_sorted_seq
noinst_PROGRAMS += sandbox/forest
noinst_PROGRAMS += sandbox/sort_pairs
noinst_PROGRAMS += sandbox/unique_pairs
#noinst_PROGRAMS += sandbox/io_into_shared_memory
noinst_PROGRAMS += tests/test_combinations
noinst_PROGRAMS += tests/test_dynamic
noinst_PROGRAMS += tests/test_read
noinst_PROGRAMS += tests/test_read_fasta
noinst_PROGRAMS += tests/test_read_serial
if HAVE_UTC_ITER
noinst_PROGRAMS += sandbox/brute_force_iter
endif

apps_align_SOURCES = apps/align.cc
apps_ccd_SOURCES = apps/ccd.c

sandbox_blast_SOURCES                      = sandbox/blast.cc
sandbox_brute_force_SOURCES                = sandbox/brute_force.cc
sandbox_brute_force_armci_SOURCES          = sandbox/brute_force_armci.cc
sandbox_brute_force_binning_SOURCES        = sandbox/brute_force_binning.cc
sandbox_brute_force_costs_SOURCES          = sandbox/brute_force.cc
sandbox_brute_force_iter_SOURCES           = sandbox/brute_force_iter.cc
sandbox_brute_force_multiple_SOURCES       = sandbox/brute_force.cc
sandbox_brute_force_pairs_SOURCES          = sandbox/brute_force_pairs.cc
sandbox_brute_force_rr_SOURCES             = sandbox/brute_force.cc
sandbox_brute_force_rr_sorted_SOURCES      = sandbox/brute_force.cc
sandbox_brute_force_rr_sorted_seq_SOURCES  = sandbox/brute_force.cc
sandbox_brute_force_serial_SOURCES         = sandbox/brute_force_serial.cc
sandbox_brute_force_sorted_SOURCES         = sandbox/brute_force.cc
sandbox_brute_force_sorted_seq_SOURCES     = sandbox/brute_force.cc
sandbox_forest_SOURCES                     = sandbox/forest.cc
sandbox_sort_pairs_SOURCES                 = sandbox/sort_pairs.cc
sandbox_unique_pairs_SOURCES               = sandbox/unique_pairs.cc
#sandbox_io_into_shared_memory_SOURCES      = sandbox/io_into_shared_memory.c
tests_test_combinations_SOURCES            = tests/test_combinations.cc
tests_test_dynamic_SOURCES                 = tests/test_dynamic.c
tests_test_read_SOURCES                    = tests/test_read.cc
tests_test_read_fasta_SOURCES              = tests/test_read_fasta.cc
tests_test_read_serial_SOURCES             = tests/test_read_serial.cc

sandbox_brute_force_costs_CPPFLAGS         = $(AM_CPPFLAGS) -DDUMP_COSTS=1
sandbox_brute_force_multiple_CPPFLAGS      = $(AM_CPPFLAGS) -DMULTIPLE_PAIRS_PER_TASK=1
sandbox_brute_force_rr_CPPFLAGS            = $(AM_CPPFLAGS) -DROUND_ROBIN_SEQUENCE_PAIRS=1
sandbox_brute_force_rr_sorted_CPPFLAGS     = $(AM_CPPFLAGS) -DROUND_ROBIN_SEQUENCE_PAIRS=1 -DSORT_TASKS_LOCALLY=1
sandbox_brute_force_rr_sorted_seq_CPPFLAGS = $(AM_CPPFLAGS) -DROUND_ROBIN_SEQUENCE_PAIRS=1 -DSORT_SEQUENCES=1
sandbox_brute_force_sorted_CPPFLAGS        = $(AM_CPPFLAGS) -DSORT_TASKS_LOCALLY=1
sandbox_brute_force_sorted_seq_CPPFLAGS    = $(AM_CPPFLAGS) -DSORT_SEQUENCES=1
