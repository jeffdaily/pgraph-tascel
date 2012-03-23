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

AM_LDFLAGS += $(TASCEL_LDFLAGS)
AM_LDFLAGS += $(GA_LDFLAGS)
AM_LDFLAGS += $(PT_MPI_LDFLAGS)

LDADD += libpgraph.la
LDADD += $(TASCEL_LIBS)
LDADD += $(GA_LIBS)
LDADD += $(PT_MPI_LIBS)
LDADD += $(GA_FLIBS)

noinst_LTLIBRARIES = libpgraph.la

libpgraph_la_SOURCES =
libpgraph_la_SOURCES += src/dynamic.c
libpgraph_la_SOURCES += src/dynamic.h

noinst_PROGRAMS =
noinst_PROGRAMS += sandbox/brute_force
noinst_PROGRAMS += sandbox/io_into_shared_memory
noinst_PROGRAMS += tests/test_dynamic
sandbox_brute_force_SOURCES           = sandbox/brute_force.cc
sandbox_io_into_shared_memory_SOURCES = sandbox/io_into_shared_memory.c
tests_test_dynamic_SOURCES            = tests/test_dynamic.c
