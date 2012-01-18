#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([HYpHMC], [0.1.0], [john.bulava@cern.ch])
AC_CONFIG_SRCDIR([lib/AnalyzerObservableGaussianWeightEstimate.C])
AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

AM_PROG_AS

AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([arpack], [main],[],AC_CHECK_LIB([arpack_Linux], [main], [], AC_MSG_ERROR([Library not Found])))
AC_CHECK_LIB([blas], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([cln], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([fftw3], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([fftw3_threads], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([gfortran], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([gomp], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([gslcblas], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([lapack], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([m], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([numa], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([pthread], [main],[],AC_MSG_ERROR([Library not Found]))
AC_CHECK_LIB([xml2], [main],[],AC_MSG_ERROR([Library not Found]))

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strrchr strstr])

AC_CONFIG_SUBDIRS(lib)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(mainprogs/Makefile)

AC_OUTPUT
