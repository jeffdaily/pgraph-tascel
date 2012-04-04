ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS =
AM_CXXFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD =

AM_CPPFLAGS += $(TASCEL_CPPFLAGS)
AM_CPPFLAGS += $(GA_CPPFLAGS)
AM_CPPFLAGS += $(PT_MPI_CPPFLAGS)
AM_CPPFLAGS += -I$(top_srcdir)/src
AM_CPPFLAGS += -I$(top_srcdir)/src/stdcomb
AM_CPPFLAGS += -I$(top_srcdir)/src/stdcomb/stdcomb

AM_LDFLAGS += $(TASCEL_LDFLAGS)
AM_LDFLAGS += $(GA_LDFLAGS)
AM_LDFLAGS += $(PT_MPI_LDFLAGS)

LDADD += libpgraph.la
LDADD += libstdcomb.la
LDADD += $(TASCEL_LIBS)
LDADD += $(GA_LIBS)
LDADD += $(PT_MPI_LIBS)
LDADD += $(GA_FLIBS)

noinst_LTLIBRARIES =
noinst_LTLIBRARIES += libpgraph.la
noinst_LTLIBRARIES += libstdcomb.la

libpgraph_la_SOURCES =
libpgraph_la_SOURCES += src/dynamic.c
libpgraph_la_SOURCES += src/dynamic.h

libstdcomb_la_SOURCES =
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigIntegerAlgorithms.cpp
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigIntegerAlgorithms.h
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigInteger.cpp
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigInteger.h
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigIntegerLibrary.h
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigIntegerUtils.cpp
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigIntegerUtils.h
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigUnsigned.cpp
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigUnsigned.h
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigUnsignedInABase.cpp
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/BigUnsignedInABase.h
libstdcomb_la_SOURCES += src/stdcomb/BigInteger/NumberlikeArray.h
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/CombFromRepSet.cpp
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/CombFromRepSet.h
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/combination.h
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/CombWithRep.cpp
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/CombWithRep.h
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/FindCombByIdx.cpp
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/FindCombByIdx.h
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/FindTotalComb.cpp
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/FindTotalComb.h
libstdcomb_la_SOURCES += src/stdcomb/stdcomb/IndexCombination.h

noinst_PROGRAMS =
noinst_PROGRAMS += sandbox/brute_force
noinst_PROGRAMS += sandbox/forest
noinst_PROGRAMS += sandbox/io_into_shared_memory
noinst_PROGRAMS += src/stdcomb/Benchmark/BenchmarkMain
noinst_PROGRAMS += src/stdcomb/CombFromRepSet/CombFromRepSetMain
noinst_PROGRAMS += src/stdcomb/CombWithRep/CombWithRepMain
noinst_PROGRAMS += src/stdcomb/FindComb/FindComb
noinst_PROGRAMS += src/stdcomb/IntComb/IntComb
noinst_PROGRAMS += tests/test_dynamic

sandbox_brute_force_SOURCES                           = sandbox/brute_force.cc
sandbox_forest_SOURCES                                = sandbox/forest.cc
sandbox_io_into_shared_memory_SOURCES                 = sandbox/io_into_shared_memory.c
src_stdcomb_Benchmark_BenchmarkMain_SOURCES           = src/stdcomb/Benchmark/BenchmarkMain.cpp
src_stdcomb_CombFromRepSet_CombFromRepSetMain_SOURCES = src/stdcomb/CombFromRepSet/CombFromRepSetMain.cpp
src_stdcomb_CombWithRep_CombWithRepMain_SOURCES       = src/stdcomb/CombWithRep/CombWithRepMain.cpp
src_stdcomb_FindComb_FindComb_SOURCES                 = src/stdcomb/FindComb/FindComb.cpp
src_stdcomb_IntComb_IntComb_SOURCES                   = src/stdcomb/IntComb/IntComb.cpp
tests_test_dynamic_SOURCES                            = tests/test_dynamic.c
