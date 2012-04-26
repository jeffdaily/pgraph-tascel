# TASCEL_REQUIRE
# --------------
AC_DEFUN([TASCEL_REQUIRE], [
TASCEL_LIBS=
TASCEL_LDFLAGS=
TASCEL_CPPFLAGS=
AC_ARG_WITH([tascel],
    [AS_HELP_STRING([--with-tascel[[=ARG]]],
        [specify location of tascel install and/or other flags])],
    [],
    [with_tascel=yes])
AS_CASE([$with_tascel],
    [yes],  [],
    [no],   [],
            [GA_ARG_PARSE(
                [with_tascel],
                [TASCEL_LIBS],
                [TASCEL_LDFLAGS],
                [TASCEL_CPPFLAGS])])
happy=yes
# Check for header.
ga_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $TASCEL_CPPFLAGS"
AC_CHECK_HEADER([tascel.h], [happy=yes], [happy=no])
CPPFLAGS="$ga_save_CPPFLAGS"
AS_IF([test "x$happy" = xno],
    [AC_MSG_FAILURE([could not locate tascel])])
# Check for library.
AC_MSG_CHECKING([for libtascel])
## add -ltascel if not already part of TASCEL_LIBS
AS_CASE([$TASCEL_LIBS],
    [*-ltascel*], [],
    [TASCEL_LIBS="$TASCEL_LIBS -ltascel"])
ga_save_LIBS="$LIBS"
ga_save_LDFLAGS="$LDFLAGS"
ga_save_CPPFLAGS="$CPPFLAGS"
LIBS="$TASCEL_LIBS $GA_LIBS $PT_MPI_LIBS $LIBS $GA_FLIBS"
LDFLAGS="$TASCEL_LDFLAGS $GA_LDFLAGS $PT_MPI_LDFLAGS $LDFLAGS"
CPPFLAGS="$TASCEL_CPPFLAGS $GA_CPPFLAGS $PT_MPI_CPPFLAGS $CPPFLAGS"
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#include <mpi.h>
#include <tascel.h>
#include <vector>
using namespace std;
using namespace tascel;
]], [[
  int dummy_argc;
  char **dummy_argv;
  ProcRank me;
  ProcRank nproc;
  int provided;
  const unsigned int numIntraRanks=1;
  MPI_Init_thread(&dummy_argc, &dummy_argv, MPI_THREAD_MULTIPLE, &provided);
  TascelConfig::initialize(numIntraRanks, MPI_COMM_WORLD);
  me = theTwoSided().getProcRank();
  nproc = theTwoSided().numProcRanks();
  TascelConfig::finalize();
  MPI_Finalize();
  return 0;
]])], [happy=yes], [happy=no])
AC_MSG_RESULT([$happy])
LIBS="$ga_save_LIBS"
LDFLAGS="$ga_save_LDFLAGS"
CPPFLAGS="$ga_save_CPPFLAGS"
AS_IF([test "x$happy" = xno],
    [AC_MSG_FAILURE([could not locate tascel])])
AC_SUBST([TASCEL_LIBS])
AC_SUBST([TASCEL_LDFLAGS])
AC_SUBST([TASCEL_CPPFLAGS])
# whether tascel has UniformTaskCollSplitIter
AC_MSG_CHECKING([whether libtascel has UniformTaskCollSplitIter])
ga_save_LIBS="$LIBS"
ga_save_LDFLAGS="$LDFLAGS"
ga_save_CPPFLAGS="$CPPFLAGS"
LIBS="$TASCEL_LIBS $GA_LIBS $PT_MPI_LIBS $LIBS $GA_FLIBS"
LDFLAGS="$TASCEL_LDFLAGS $GA_LDFLAGS $PT_MPI_LDFLAGS $LDFLAGS"
CPPFLAGS="$TASCEL_CPPFLAGS $GA_CPPFLAGS $PT_MPI_CPPFLAGS $CPPFLAGS"
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollIter.h>
#include <vector>
using namespace std;
using namespace tascel;
static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) { }
typedef int task_description;
]], [[
  int dummy_argc;
  char **dummy_argv;
  ProcRank me;
  ProcRank nproc;
  int provided;
  const unsigned int numIntraRanks=1;
  MPI_Init_thread(&dummy_argc, &dummy_argv, MPI_THREAD_MULTIPLE, &provided);
  TascelConfig::initialize(numIntraRanks, MPI_COMM_WORLD);
  me = theTwoSided().getProcRank();
  nproc = theTwoSided().numProcRanks();
  TslFuncRegTbl *frt = new TslFuncRegTbl();
  TslFunc tf = frt->add(alignment_task);
  TaskCollProps props;
  props.functions(tf, frt)
      .taskSize(sizeof(task_description));
  UniformTaskCollIter utc(props, NULL, 0, -1, 0);
  TascelConfig::finalize();
  MPI_Finalize();
  return 0;
]])], [happy=yes], [happy=no])
AC_MSG_RESULT([$happy])
AM_CONDITIONAL([HAVE_UTC_ITER], [test "x$happy" = xyes])
])dnl
