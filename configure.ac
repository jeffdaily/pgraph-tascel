#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([parallel sequence graph construction using work stealing],
        [0.1], [jeff.daily@pnnl.gov], pgraph-tascel)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([dynamic.h])

AC_LANG([C++])

# Checks for programs.
AC_PROG_CXX([mpic++ mpicxx mpiCC sxmpic++ hcp mpxlC_r mpxlC mpixlcxx_r mpixlcxx mpg++ mpc++ mpCC cmpic++ mpiFCC CC])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

LT_INIT([disable-shared])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
