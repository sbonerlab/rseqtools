#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([rseqtools], [0.7.0], [david.z.chen@outlook.com])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])

#------------------------------------------------------------------------------
# Initialize automake
#------------------------------------------------------------------------------
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------
AC_PROG_CC

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------
AC_CHECK_LIB([m], [log], [], [AC_MSG_ERROR([Cannot find standard math library])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Cannot find cblas library])])
AC_CHECK_LIB([gsl], [gsl_ran_hypergeometric_pdf], [], [AC_MSG_ERROR([Cannot find gsl library])])
AC_CHECK_LIB([bios], [needMem], [], [AC_MSG_ERROR([Cannot find bios library])])
AC_CHECK_LIB([mrf], [mrf_init], [], [AC_MSG_ERROR([Cannot find mrf library])])

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------
AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

#------------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------------
AC_HEADER_STDBOOL

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------
AC_CHECK_FUNCS([strchr strdup strstr])

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
