# GA_SEARCH_LIBS(function, search-LIBS,
#                [action-if-found], [action-if-not-found],
#                [other-libraries], [preamble])
# --------------------------------------------------------
# Like AC_SEARCH_LIBS, but allow an optional preamble.
AC_DEFUN([GA_SEARCH_LIBS], [
AS_VAR_PUSHDEF([ga_search], [ga_cv_search_$1])dnl
AC_CACHE_CHECK([for library containing $1], [ga_search],
    [ga_func_search_save_LIBS=$LIBS
     AC_LANG_CONFTEST([AC_LANG_CALL([$6], [$1])])
     for ga_lib in '' $2
     do
        if test -z "$ga_lib"; then
            ga_res="none required"
        else
            ga_res=-l$ga_lib
            LIBS="-l$ga_lib $5 $ga_func_search_save_LIBS"
        fi
        AC_LINK_IFELSE([], [AS_VAR_SET([ga_search], [$ga_res])])
        AS_VAR_SET_IF([ga_search], [break])
     done
     AS_VAR_SET_IF([ga_search], , [AS_VAR_SET([ga_search], [no])])
     rm conftest.$ac_ext
     LIBS=$ga_func_search_save_LIBS])
AS_VAR_COPY([ga_res], [ga_search])
AS_IF([test "$ga_res" != no],
    [test "$ga_res" = "none required" || LIBS="$ga_res $LIBS"
     $3],
    [$4])
AS_VAR_POPDEF([ga_search])dnl
])


# GA_CHECK_PACKAGE(pkg, header, library, function, [extra-libs],
#                  [action-if-found], [action-if-not-found],
#                  [function-preamble])
# --------------------------------------------------------------
#
AC_DEFUN([GA_CHECK_PACKAGE], [
AS_VAR_PUSHDEF([HAVE_PKG],    m4_toupper(m4_translit([HAVE_$1], [-.], [__])))
AS_VAR_PUSHDEF([PKG_LIBS],    m4_toupper(m4_translit([$1_LIBS], [-.], [__])))
AS_VAR_PUSHDEF([PKG_LDFLAGS], m4_toupper(m4_translit([$1_LDFLAGS], [-.], [__])))
AS_VAR_PUSHDEF([PKG_CPPFLAGS],m4_toupper(m4_translit([$1_CPPFLAGS], [-.], [__])))
AS_VAR_SET([PKG_LIBS],[])
AS_VAR_SET([PKG_LDFLAGS],[])
AS_VAR_SET([PKG_CPPFLAGS],[])
AC_ARG_WITH([$1],
    [AS_HELP_STRING([--with-$1[[=ARG]]],
        [specify location of $1 install and/or other flags])],
    [],
    [with_$1=yes])
AS_CASE([$with_$1],
    [yes],  [],
    [no],   [],
            [GA_ARG_PARSE(
                [with_$1],
                [PKG_LIBS],
                [PKG_LDFLAGS],
                [PKG_CPPFLAGS])])
# Check for header.
ga_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $PKG_CPPFLAGS"
AC_CHECK_HEADER([$2], [], [$7])
CPPFLAGS="$ga_save_CPPFLAGS"
# Check for library.
ga_save_LIBS="$LIBS"; LIBS="$PKG_LIBS $LIBS"
ga_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LDFLAGS $PKG_LDFLAGS"
ga_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $PKG_CPPFLAGS"
GA_SEARCH_LIBS([$4], [$3], [], [], [$5], [$8])
LIBS="$ga_save_LIBS"
LDFLAGS="$ga_save_LDFLAGS"
CPPFLAGS="$ga_save_CPPFLAGS"
AS_IF([test "x$ga_cv_search_$4" != xno],
    [$6
     AS_IF([test "x$ga_cv_search_$4" != "xnone required"],
        [AS_VAR_APPEND([PKG_LIBS], [$ga_cv_search_$4])])
     AC_DEFINE([HAVE_PKG], [1], [set to 1 if we have the indicated package])
     AC_SUBST(PKG_LIBS)
     AC_SUBST(PKG_LDFLAGS)
     AC_SUBST(PKG_CPPFLAGS)],
    [$7])
AM_CONDITIONAL(HAVE_PKG, [test "x$ga_cv_search_$4" != xno])
AS_VAR_POPDEF([HAVE_PKG])
AS_VAR_POPDEF([PKG_LIBS])
AS_VAR_POPDEF([PKG_LDFLAGS])
AS_VAR_POPDEF([PKG_CPPFLAGS])
])dnl
