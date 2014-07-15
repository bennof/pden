
AC_PREREQ([2.68])
AC_INIT([PDEN], [0.3], [pden@benjaminfalkner.de])
AM_INIT_AUTOMAKE([0.3 foreign])
AC_CONFIG_SRCDIR([src/pden.h])
AC_CONFIG_HEADERS([src/config.h])

ACLOCAL="$ACLOCAL $ACLOCAL_FLAGS"


# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB


# Threading 
AC_OPENMP

# Checks for libraries.
AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c_3d])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset sqrt strstr])

# Debug
AC_ARG_ENABLE(debug,AC_HELP_STRING(--enable-debug, create a debug build),found_debug=\$enable_debug)
AM_CONDITIONAL(IS_DEBUG, test "x\$found_debug" = xyes)

# Precision
AC_ARG_ENABLE(debug,AC_HELP_STRING(--enable-single, create a single precision build),found_single=\$enable_single)
AM_CONDITIONAL(IS_SINGLE, test "x\$found_single" = xyes)

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
