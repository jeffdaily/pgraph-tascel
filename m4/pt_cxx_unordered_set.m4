# PT_CXX_UNORDERED_SET
# --------------------
# Look for <unordered_set> then <tr1/unordered_set>. Defines
# HAVE_CXX_UNORDERED_SET and HAVE_CXX_TR1_UNORDERED_SET, respectively.
AC_DEFUN([PT_CXX_UNORDERED_SET], [
AC_CACHE_CHECK([for unordered_set],
    [pt_cv_cxx_unordered_set],
    [AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM(
[[#include <unordered_set>
using std::unordered_set;]],
[[unordered_set<int> my_test_set;]])],
        [pt_cv_cxx_unordered_set=yes; pt_cv_cxx_unordered_set_val=1],
        [pt_cv_cxx_unordered_set=no; pt_cv_cxx_unordered_set_val=0])])

AC_DEFINE_UNQUOTED([HAVE_CXX_UNORDERED_SET],
                   [$pt_cv_cxx_unordered_set_val],
                   [Define if unordered_set is present.])

AC_CACHE_CHECK([for tr1/unordered_set],
    [pt_cv_cxx_tr1_unordered_set],
    [AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM(
[[#include <tr1/unordered_set>
using std::tr1::unordered_set;]],
[[unordered_set<int> my_test_set;]])],
        [pt_cv_cxx_tr1_unordered_set=yes; pt_cv_cxx_tr1_unordered_set_val=1],
        [pt_cv_cxx_tr1_unordered_set=no; pt_cv_cxx_tr1_unordered_set_val=0])])

AC_DEFINE_UNQUOTED([HAVE_CXX_TR1_UNORDERED_SET],
                   [$pt_cv_cxx_tr1_unordered_set_val],
                   [Define if tr1/unordered_set is present.])
])
