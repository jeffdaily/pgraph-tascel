ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS =
AM_CXXFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD =

AM_CPPFLAGS += $(GA_CPPFLAGS)
AM_CPPFLAGS += $(PT_MPI_CPPFLAGS)
AM_CPPFLAGS += -I$(top_srcdir)/src

AM_LDFLAGS += $(GA_LDFLAGS)
AM_LDFLAGS += $(PT_MPI_LDFLAGS)

LDADD += libpgraph.la
LDADD += $(GA_LIBS)
LDADD += $(PT_MPI_LIBS)
LDADD += $(GA_FLIBS)

lib_LTLIBRARIES = libpgraph.la

libpgraph_la_SOURCES =
libpgraph_la_SOURCES += src/dynamic.c
libpgraph_la_SOURCES += src/dynamic.h

check_PROGRAMS = test_dynamic
test_dynamic_SOURCES = test_dynamic.c
