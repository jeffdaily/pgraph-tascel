# PT_WITH_MPI
# -----------
# Establishes all things related to the Message Passing Interface.
# This includes the compilers to use (either standard or MPI wrappers)
# or the proper linker flags (-L), libs (-l) or preprocessor directives (-I).
# Yes, it's a beefy AC macro, but because when MPI is desired it replaces the
# usual compiler the order here is necessary and it is all interdependent.
AC_DEFUN([PT_WITH_MPI], [
# PT_MPI_* vars might exist in environment, but they are really internal.
# Reset them.
PT_MPI_LIBS=
PT_MPI_LDFLAGS=
PT_MPI_CPPFLAGS=
AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi[[=ARG]]],
        [select MPI as the messaging library (default); leave ARG blank to use MPI compiler wrappers])],
    [],
    [with_mpi=yes])
with_mpi_need_parse=no
AS_CASE([$with_mpi],
    [yes],  [with_mpi_wrappers=yes],
    [no],   [],
    [*],    [with_mpi_need_parse=yes])
dnl postpone parsing with_mpi until we know sizeof(void*)
dnl AS_IF([test x$with_mpi_need_parse = xyes],
dnl     [PT_ARG_PARSE([with_mpi], [PT_MPI_LIBS], [PT_MPI_LDFLAGS], [PT_MPI_CPPFLAGS])])
AC_SUBST([PT_MPI_LIBS])
AC_SUBST([PT_MPI_LDFLAGS])
AC_SUBST([PT_MPI_CPPFLAGS])
])dnl
