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
if HAVE_SSE2
AM_CPPFLAGS += -I$(top_srcdir)/contrib/Complete-Striped-Smith-Waterman-Library/src
endif

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
libpgraph_la_SOURCES += src/alignment.cpp
libpgraph_la_SOURCES += src/alignment.hpp
libpgraph_la_SOURCES += src/AlignStats.hpp
libpgraph_la_SOURCES += src/blosum/blosum40.h
libpgraph_la_SOURCES += src/blosum/blosum45.h
libpgraph_la_SOURCES += src/blosum/blosum62.h
libpgraph_la_SOURCES += src/blosum/blosum75.h
libpgraph_la_SOURCES += src/blosum/blosum80.h
libpgraph_la_SOURCES += src/blosum/blosum90.h
libpgraph_la_SOURCES += src/bucket.cpp
libpgraph_la_SOURCES += src/bucket.hpp
libpgraph_la_SOURCES += src/combinations.c
libpgraph_la_SOURCES += src/combinations.h
libpgraph_la_SOURCES += src/constants.h
libpgraph_la_SOURCES += src/csequence.c
libpgraph_la_SOURCES += src/csequence.h
libpgraph_la_SOURCES += src/io.cpp
libpgraph_la_SOURCES += src/io.hpp
libpgraph_la_SOURCES += src/mpix.cpp
libpgraph_la_SOURCES += src/mpix.hpp
libpgraph_la_SOURCES += src/Parameters.cpp
libpgraph_la_SOURCES += src/Parameters.hpp
libpgraph_la_SOURCES += src/param.cpp
libpgraph_la_SOURCES += src/param.hpp
libpgraph_la_SOURCES += src/pthread_fixes.h
libpgraph_la_SOURCES += src/Sequence.cpp
libpgraph_la_SOURCES += src/Sequence.hpp
libpgraph_la_SOURCES += src/SequenceDatabaseArmci.cpp
libpgraph_la_SOURCES += src/SequenceDatabaseArmci.hpp
libpgraph_la_SOURCES += src/SequenceDatabaseGArray.cpp
libpgraph_la_SOURCES += src/SequenceDatabaseGArray.hpp
libpgraph_la_SOURCES += src/SequenceDatabase.cpp
libpgraph_la_SOURCES += src/SequenceDatabase.hpp
libpgraph_la_SOURCES += src/stree.cpp
libpgraph_la_SOURCES += src/stree.hpp
libpgraph_la_SOURCES += src/SuffixBuckets.cpp
libpgraph_la_SOURCES += src/SuffixBuckets.hpp
libpgraph_la_SOURCES += src/SuffixTree.cpp
libpgraph_la_SOURCES += src/SuffixTree.hpp
libpgraph_la_SOURCES += src/TreeStats.hpp
libpgraph_la_SOURCES += src/timer.h
if HAVE_SSE2
libpgraph_la_SOURCES += contrib/Complete-Striped-Smith-Waterman-Library/src/ssw.c
libpgraph_la_SOURCES += contrib/Complete-Striped-Smith-Waterman-Library/src/ssw.h
libpgraph_la_SOURCES += contrib/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.cpp
libpgraph_la_SOURCES += contrib/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h
endif

bin_PROGRAMS =
bin_PROGRAMS += apps/align
bin_PROGRAMS += apps/align_counter
bin_PROGRAMS += apps/align_counter_length_filter
bin_PROGRAMS += apps/align_garray
bin_PROGRAMS += apps/align_iter
bin_PROGRAMS += apps/align_iter_length_filter
bin_PROGRAMS += apps/align_length_filter
bin_PROGRAMS += apps/align_length_filter_noalign
bin_PROGRAMS += apps/align_tree
bin_PROGRAMS += apps/align_tree_noalign
bin_PROGRAMS += apps/align_tree_dumb_noalign
bin_PROGRAMS += apps/align_tree_dup
bin_PROGRAMS += apps/align_tree_dynamic
bin_PROGRAMS += apps/align_tree_dynamic_dup
bin_PROGRAMS += apps/align_tree_dynamic_callback
bin_PROGRAMS += apps/align_tree_dynamic_callback_dup
bin_PROGRAMS += apps/align_tree_hybrid
bin_PROGRAMS += apps/align_tree_hybrid_noalign
bin_PROGRAMS += apps/align_tree_hybrid_sorted
bin_PROGRAMS += apps/align_tree_hybrid_sorted_noalign
bin_PROGRAMS += apps/align_tree_hybrid_dup
bin_PROGRAMS += apps/align_tree_hybrid_callback
bin_PROGRAMS += apps/align_tree_hybrid_callback_sorted
bin_PROGRAMS += apps/align_tree_hybrid_callback_sorted_noalign
bin_PROGRAMS += apps/align_tree_hybrid_callback_dup
bin_PROGRAMS += apps/ccd
bin_PROGRAMS += apps/ccd2
if HAVE_SSE2
bin_PROGRAMS += apps/align_ssw
bin_PROGRAMS += apps/align_counter_ssw
bin_PROGRAMS += apps/align_iter_ssw
bin_PROGRAMS += apps/align_tree_ssw
if HAVE_ZLIB
bin_PROGRAMS += apps/ssw_test
endif
endif

noinst_PROGRAMS =
#noinst_PROGRAMS += sandbox/blast
#noinst_PROGRAMS += sandbox/brute_force
#if HAVE_ARMCI
#noinst_PROGRAMS += sandbox/brute_force_armci
#endif
#noinst_PROGRAMS += sandbox/brute_force_binning
#noinst_PROGRAMS += sandbox/brute_force_costs
#noinst_PROGRAMS += sandbox/brute_force_multiple
#noinst_PROGRAMS += sandbox/brute_force_pairs
#noinst_PROGRAMS += sandbox/brute_force_rr
#noinst_PROGRAMS += sandbox/brute_force_rr_sorted
#noinst_PROGRAMS += sandbox/brute_force_rr_sorted_seq
#noinst_PROGRAMS += sandbox/brute_force_serial
#noinst_PROGRAMS += sandbox/brute_force_sorted
#noinst_PROGRAMS += sandbox/brute_force_sorted_seq
#noinst_PROGRAMS += sandbox/forest
#noinst_PROGRAMS += sandbox/sort_pairs
#noinst_PROGRAMS += sandbox/unique_pairs
##noinst_PROGRAMS += sandbox/io_into_shared_memory
noinst_PROGRAMS += tests/test_align
noinst_PROGRAMS += tests/test_combinations
noinst_PROGRAMS += tests/test_dynamic
noinst_PROGRAMS += tests/test_kcombo2
noinst_PROGRAMS += tests/test_read
noinst_PROGRAMS += tests/test_read_fasta
noinst_PROGRAMS += tests/test_read_serial
if HAVE_SSE2
noinst_PROGRAMS += tests/test_ssw
endif
#if HAVE_UTC_ITER
#noinst_PROGRAMS += sandbox/brute_force_iter
#endif

apps_align_SOURCES = apps/align.cpp
apps_align_counter_SOURCES = apps/align_counter.cpp
apps_align_counter_length_filter_SOURCES = apps/align_counter.cpp
apps_align_garray_SOURCES = apps/align.cpp
apps_align_iter_SOURCES = apps/align.cpp
apps_align_iter_length_filter_SOURCES = apps/align.cpp
apps_align_length_filter_SOURCES = apps/align.cpp
apps_align_length_filter_noalign_SOURCES = apps/align.cpp
apps_align_tree_SOURCES = apps/align_tree.cpp
apps_align_tree_noalign_SOURCES = apps/align_tree.cpp
apps_align_tree_dumb_noalign_SOURCES = apps/align_tree.cpp
apps_align_tree_dup_SOURCES = apps/align_tree.cpp
apps_align_tree_dynamic_SOURCES = apps/align_tree_dynamic.cpp
apps_align_tree_dynamic_dup_SOURCES = apps/align_tree_dynamic.cpp
apps_align_tree_dynamic_callback_SOURCES = apps/align_tree_dynamic.cpp
apps_align_tree_dynamic_callback_dup_SOURCES = apps/align_tree_dynamic.cpp
apps_align_tree_hybrid_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_noalign_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_sorted_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_sorted_noalign_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_dup_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_callback_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_callback_sorted_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_callback_sorted_noalign_SOURCES = apps/align_tree_hybrid.cpp
apps_align_tree_hybrid_callback_dup_SOURCES = apps/align_tree_hybrid.cpp
apps_ccd_SOURCES = apps/ccd.cpp
apps_ccd2_SOURCES = apps/ccd2.cpp
apps_align_ssw_SOURCES = apps/align.cpp
apps_align_counter_ssw_SOURCES = apps/align_counter.cpp
apps_align_iter_ssw_SOURCES = apps/align.cpp
apps_align_tree_ssw_SOURCES = apps/align_tree.cpp
apps_ssw_test_SOURCES =
apps_ssw_test_SOURCES += contrib/Complete-Striped-Smith-Waterman-Library/src/main.c
apps_ssw_test_SOURCES += contrib/Complete-Striped-Smith-Waterman-Library/src/kseq.h

apps_align_iter_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_ITER
apps_align_iter_length_filter_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_ITER -DLENGTH_FILTER
apps_align_length_filter_CPPFLAGS = $(AM_CPPFLAGS) -DLENGTH_FILTER
apps_align_length_filter_noalign_CPPFLAGS = $(AM_CPPFLAGS) -DLENGTH_FILTER -DNOALIGN
apps_align_garray_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_GARRAY
apps_align_ssw_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_SSW
apps_align_counter_length_filter_CPPFLAGS = $(AM_CPPFLAGS) -DLENGTH_FILTER
apps_align_counter_ssw_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_SSW
apps_align_iter_ssw_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_SSW
apps_align_tree_noalign_CPPFLAGS = $(AM_CPPFLAGS) -DNOALIGN
apps_align_tree_dumb_noalign_CPPFLAGS = $(AM_CPPFLAGS) -DNOALIGN -DDUMB
apps_align_tree_ssw_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_SSW
apps_align_tree_dup_CPPFLAGS = $(AM_CPPFLAGS) -DGLOBAL_DUPLICATES
apps_align_tree_dynamic_dup_CPPFLAGS = $(AM_CPPFLAGS) -DGLOBAL_DUPLICATES
apps_align_tree_dynamic_callback_CPPFLAGS = $(AM_CPPFLAGS) -DCALLBACK
apps_align_tree_dynamic_callback_dup_CPPFLAGS = $(AM_CPPFLAGS) -DCALLBACK -DGLOBAL_DUPLICATES
apps_align_tree_hybrid_noalign_CPPFLAGS = $(AM_CPPFLAGS) -DNOALIGN
apps_align_tree_hybrid_sorted_CPPFLAGS = $(AM_CPPFLAGS) -DSORT_BUCKETS
apps_align_tree_hybrid_sorted_noalign_CPPFLAGS = $(AM_CPPFLAGS) -DNOALIGN -DSORT_BUCKETS
apps_align_tree_hybrid_dup_CPPFLAGS = $(AM_CPPFLAGS) -DGLOBAL_DUPLICATES
apps_align_tree_hybrid_callback_CPPFLAGS = $(AM_CPPFLAGS) -DCALLBACK
apps_align_tree_hybrid_callback_sorted_CPPFLAGS = $(AM_CPPFLAGS) -DCALLBACK -DSORT_BUCKETS
apps_align_tree_hybrid_callback_sorted_noalign_CPPFLAGS = $(AM_CPPFLAGS) -DNOALIGN -DCALLBACK -DSORT_BUCKETS
apps_align_tree_hybrid_callback_dup_CPPFLAGS = $(AM_CPPFLAGS) -DCALLBACK -DGLOBAL_DUPLICATES

apps_ssw_test_LDFLAGS = $(ZLIB_LDFLAGS) $(AM_LDFLAGS)

apps_ssw_test_LDADD = $(LDADD) $(ZLIB_LIBS) $(LIBS)

#sandbox_blast_SOURCES                      = sandbox/blast.cpp
#sandbox_brute_force_SOURCES                = sandbox/brute_force.cpp
#sandbox_brute_force_armci_SOURCES          = sandbox/brute_force_armci.cpp
#sandbox_brute_force_binning_SOURCES        = sandbox/brute_force_binning.cpp
#sandbox_brute_force_costs_SOURCES          = sandbox/brute_force.cpp
#sandbox_brute_force_iter_SOURCES           = sandbox/brute_force_iter.cpp
#sandbox_brute_force_multiple_SOURCES       = sandbox/brute_force.cpp
#sandbox_brute_force_pairs_SOURCES          = sandbox/brute_force_pairs.cpp
#sandbox_brute_force_rr_SOURCES             = sandbox/brute_force.cpp
#sandbox_brute_force_rr_sorted_SOURCES      = sandbox/brute_force.cpp
#sandbox_brute_force_rr_sorted_seq_SOURCES  = sandbox/brute_force.cpp
#sandbox_brute_force_serial_SOURCES         = sandbox/brute_force_serial.cpp
#sandbox_brute_force_sorted_SOURCES         = sandbox/brute_force.cpp
#sandbox_brute_force_sorted_seq_SOURCES     = sandbox/brute_force.cpp
#sandbox_forest_SOURCES                     = sandbox/forest.cpp
#sandbox_sort_pairs_SOURCES                 = sandbox/sort_pairs.cpp
#sandbox_unique_pairs_SOURCES               = sandbox/unique_pairs.cpp
##sandbox_io_into_shared_memory_SOURCES      = sandbox/io_into_shared_memory.c
tests_test_align_SOURCES                   = tests/test_align.cpp
tests_test_combinations_SOURCES            = tests/test_combinations.cpp
tests_test_dynamic_SOURCES                 = tests/test_dynamic.cpp
tests_test_kcombo2_SOURCES                 = tests/test_kcombo2.cpp
tests_test_read_SOURCES                    = tests/test_read.cpp
tests_test_read_fasta_SOURCES              = tests/test_read_fasta.cpp
tests_test_read_serial_SOURCES             = tests/test_read_serial.cpp
tests_test_ssw_SOURCES                     = tests/test_ssw.cpp
#
#sandbox_brute_force_costs_CPPFLAGS         = $(AM_CPPFLAGS) -DDUMP_COSTS=1
#sandbox_brute_force_multiple_CPPFLAGS      = $(AM_CPPFLAGS) -DMULTIPLE_PAIRS_PER_TASK=1
#sandbox_brute_force_rr_CPPFLAGS            = $(AM_CPPFLAGS) -DROUND_ROBIN_SEQUENCE_PAIRS=1
#sandbox_brute_force_rr_sorted_CPPFLAGS     = $(AM_CPPFLAGS) -DROUND_ROBIN_SEQUENCE_PAIRS=1 -DSORT_TASKS_LOCALLY=1
#sandbox_brute_force_rr_sorted_seq_CPPFLAGS = $(AM_CPPFLAGS) -DROUND_ROBIN_SEQUENCE_PAIRS=1 -DSORT_SEQUENCES=1
#sandbox_brute_force_sorted_CPPFLAGS        = $(AM_CPPFLAGS) -DSORT_TASKS_LOCALLY=1
#sandbox_brute_force_sorted_seq_CPPFLAGS    = $(AM_CPPFLAGS) -DSORT_SEQUENCES=1
