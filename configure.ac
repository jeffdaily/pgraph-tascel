#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
###############################################################################
# Init autoconf
###############################################################################
AC_PREREQ([2.68])
AC_INIT([parallel sequence graph construction using work stealing],
        [0.1], [jeff.daily@pnnl.gov], pgraph-tascel)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/Sequence.cpp])

###############################################################################
# Init automake
###############################################################################
AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

###############################################################################
# Misc initial setup
###############################################################################
# MPI is required; user can specify to use MPI compilers or standard compilers
PT_WITH_MPI

###############################################################################
# C compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C compiler])
AC_MSG_NOTICE

PT_PROG_MPICC
AC_USE_SYSTEM_EXTENSIONS

AS_IF([test x$with_mpi_wrappers = xyes],
    [PT_MPI_UNWRAP],
    [GA_ARG_PARSE([with_mpi], [PT_MPI_LIBS], [PT_MPI_LDFLAGS], [PT_MPI_CPPFLAGS])])

# Checks for libraries.
#GA_REQUIRE

# Checks for header files.
AC_CHECK_HEADERS([emmintrin.h])
AC_CHECK_HEADERS([limits.h])
AC_CHECK_HEADERS([pthread.h])
AC_CHECK_HEADERS([sched.h])
AC_CHECK_HEADERS([stdlib.h])

AM_CONDITIONAL([HAVE_SSE2], [test "x$ac_cv_header_emmintrin_h" = xyes])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T

AC_CHECK_TYPES([pthread_barrier_t], [], [], [[#include <pthread.h>]])
AM_CONDITIONAL([HAVE_PTHREAD_BARRIER_T],
               [test "x$ac_cv_type_pthread_barrier_t" = xyes])
AH_BOTTOM([#if !HAVE_PTHREAD_BARRIER_T
#include "pthread_fixes.h"
#endif])

AC_CHECK_FUNCS([sched_setaffinity])
AC_CHECK_FUNCS([pthread_setaffinity_np])
AS_CASE([$ac_cv_func_sched_setaffinity:$ac_cv_func_pthread_setaffinity_np],
        [yes:yes],
        [AC_DEFINE([SET_AFFINITY], [1], [whether affinity is supported])])
dnl         [set_affinity=1],
dnl         [set_affinity=0])
dnl AC_DEFINE_UNQUOTED([SET_AFFINITY], [$set_affinity],
dnl                    [whether affinity is supported])

# Checks for library functions.
AC_SEARCH_LIBS([abs], [m])

###############################################################################
# C++ compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C++ compiler])
AC_MSG_NOTICE
AC_LANG_PUSH([C++])

PT_PROG_MPICXX

AS_IF([test x$with_mpi_wrappers = xyes], [PT_MPI_UNWRAP])

# Checks for libraries.
TASCEL_REQUIRE
GA_CHECK_PACKAGE([gmp], [gmp.h], [gmp], [__gmpz_bin_uiui])
GA_CHECK_PACKAGE([armci], [armci.h], [armci], [ARMCI_Init])
GA_CHECK_PACKAGE([zlib], [zlib.h], [z], [gzread])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_LANG_POP([C++])

###############################################################################
# Libtool setup -- no compiler/linker tests after this
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Libtool setup])
AC_MSG_NOTICE
# temporarily restore unwrapped compilers
# this works around a bug where libtool sadly relies on matching compiler
# names in order to determine features (Fortran only, I think)
# libtool doesn't recognize MPI compiler names, nor should it
AS_IF([test x$with_mpi_wrappers = xyes], [PT_MPI_UNWRAP_PUSH])
LT_INIT([disable-shared])
# and now that that's over, put the MPI compilers back
# also, the above hack incorrectly sets the base compiler as the linker
AS_IF([test x$with_mpi_wrappers = xyes],
    [PT_MPI_UNWRAP_POP
     compiler="$CC"
     LTCC="$CC"
     lt_save_CC="$CC"
     compiler_DEFAULT="$CC"
     compiler_CXX="$CXX"
     compiler_F77="$F77"])

###############################################################################
# The End
###############################################################################
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# Report on what we found.
AC_MSG_NOTICE([])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([ $PACKAGE_NAME])
AC_MSG_NOTICE([ configured as follows:])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([           PT_MPI_LIBS=$PT_MPI_LIBS])
AC_MSG_NOTICE([        PT_MPI_LDFLAGS=$PT_MPI_LDFLAGS])
AC_MSG_NOTICE([       PT_MPI_CPPFLAGS=$PT_MPI_CPPFLAGS])
AC_MSG_NOTICE([                    CC=$CC])
AS_IF([test "x$with_mpi_wrappers" = xyes], [
AC_MSG_NOTICE([          unwrapped CC=$pt_cv_mpic_naked])
])
AC_MSG_NOTICE([                CFLAGS=$CFLAGS])
AC_MSG_NOTICE([                   CXX=$CXX])
AS_IF([test "x$with_mpi_wrappers" = xyes], [
AC_MSG_NOTICE([         unwrapped CXX=$pt_cv_mpicxx_naked])
])
AC_MSG_NOTICE([              CXXFLAGS=$CXXFLAGS])
AC_MSG_NOTICE([             PT_CXXOPT=$PT_CXXOPT])
AC_MSG_NOTICE([                   CPP=$CPP])
AC_MSG_NOTICE([              CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([               LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([                  LIBS=$LIBS])
AC_MSG_NOTICE([        ARMCI_CPPFLAGS=$ARMCI_CPPFLAGS])
AC_MSG_NOTICE([         ARMCI_LDFLAGS=$ARMCI_LDFLAGS])
AC_MSG_NOTICE([            ARMCI_LIBS=$ARMCI_LIBS])
AC_MSG_NOTICE([           ARMCI_FLIBS=$ARMCI_FLIBS])
dnl AC_MSG_NOTICE([           GA_CPPFLAGS=$GA_CPPFLAGS])
dnl AC_MSG_NOTICE([            GA_LDFLAGS=$GA_LDFLAGS])
dnl AC_MSG_NOTICE([               GA_LIBS=$GA_LIBS])
dnl AC_MSG_NOTICE([              GA_FLIBS=$GA_FLIBS])
AC_MSG_NOTICE([       TASCEL_CPPFLAGS=$TASCEL_CPPFLAGS])
AC_MSG_NOTICE([        TASCEL_LDFLAGS=$TASCEL_LDFLAGS])
AC_MSG_NOTICE([           TASCEL_LIBS=$TASCEL_LIBS])
