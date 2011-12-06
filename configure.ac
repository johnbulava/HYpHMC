#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([HYpHMC], [0.1.0], [john.bulava@cern.ch])
AC_CONFIG_SRCDIR([AnalyzerObservableGaussianWeightEstimate.C])
AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-larpack':
AC_CHECK_LIB([arpack], [main])
# FIXME: Replace `main' with a function in `-lblas':
AC_CHECK_LIB([blas], [main])
# FIXME: Replace `main' with a function in `-lcblas':
AC_CHECK_LIB([cblas], [main])
# FIXME: Replace `main' with a function in `-lcln':
AC_CHECK_LIB([cln], [main])
# FIXME: Replace `main' with a function in `-lfftw3':
AC_CHECK_LIB([fftw3], [main])
# FIXME: Replace `main' with a function in `-lfftw3_threads':
AC_CHECK_LIB([fftw3_threads], [main])
# FIXME: Replace `main' with a function in `-lg2c':
AC_CHECK_LIB([g2c], [main])
# FIXME: Replace `main' with a function in `-lgfortran':
AC_CHECK_LIB([gfortran], [main])
# FIXME: Replace `main' with a function in `-lgomp':
AC_CHECK_LIB([gomp], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])
# FIXME: Replace `main' with a function in `-limf':
AC_CHECK_LIB([imf], [main])
# FIXME: Replace `main' with a function in `-liomp5':
AC_CHECK_LIB([iomp5], [main])
# FIXME: Replace `main' with a function in `-llapack':
AC_CHECK_LIB([lapack], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lmkl':
AC_CHECK_LIB([mkl], [main])
# FIXME: Replace `main' with a function in `-lmkl_core':
AC_CHECK_LIB([mkl_core], [main])
# FIXME: Replace `main' with a function in `-lmkl_intel_lp64':
AC_CHECK_LIB([mkl_intel_lp64], [main])
# FIXME: Replace `main' with a function in `-lmkl_lapack95_lp64':
AC_CHECK_LIB([mkl_lapack95_lp64], [main])
# FIXME: Replace `main' with a function in `-lmkl_sequential':
AC_CHECK_LIB([mkl_sequential], [main])
# FIXME: Replace `main' with a function in `-lnuma':
AC_CHECK_LIB([numa], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])
# FIXME: Replace `main' with a function in `-lxSSE64':
AC_CHECK_LIB([xSSE64], [main])
# FIXME: Replace `main' with a function in `-lxml2':
AC_CHECK_LIB([xml2], [main])

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strrchr strstr])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
