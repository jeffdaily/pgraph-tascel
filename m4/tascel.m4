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
ga_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $TASCEL_CPPFLAGS"
AC_CHECK_HEADER([UniformTaskCollection.h], [happy=yes], [happy=no])
CPPFLAGS="$ga_save_CPPFLAGS"
AS_IF([test "x$happy" = xno],
    [AC_MSG_FAILURE([could not locate tascel])])
# Check for library.
AC_MSG_CHECKING([for libtascel])
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
#include <ga.h>
#include <UniformTaskCollectionShared.h>
using namespace std;
using namespace tascel;
#define MAX_TASKS 100

static int me;

struct TskDscr {
  int v;
  TskDscr(int _v): v(_v) {}
};

void fn(UniformTaskCollection *coll, void *desc, int dscr_len, void *pldata,
        int pldata_len, std::vector<void*> data_bufs) {
  TskDscr *dsc = (TskDscr*)desc;
  if(dsc->v>0) {
    TskDscr dsc2(dsc->v-1);
    coll->addTask(&dsc2, sizeof(dsc2));
  }
}

void the_test() {
  TslFuncRegTbl frt;
  TslFunc tf = frt.add(fn);
  TaskCollProps props;
  props.functions(tf,frt).taskSize(sizeof(TskDscr)).maxTasks(MAX_TASKS);
  UniformTaskCollectionShared utc(props);
  if(me==0) {
    TskDscr dsc(10);
    utc.addTask(&dsc, sizeof(dsc));
  }
  utc.process();
}
]], [[
  int dummy_argc;
  char **dummy_argv;
  MPI_Init(&dummy_argc, &dummy_argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  GA_Initialize();
  the_test();
  GA_Terminate();
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
])dnl
