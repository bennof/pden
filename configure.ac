
AC_PREREQ([2.68])
AC_INIT([PDEN], [0.3], [pden@benjaminfalkner.de])
AM_INIT_AUTOMAKE([0.3 foreign])
AC_CONFIG_SRCDIR([src/pden.in.h])
AC_CONFIG_HEADERS([src/config.h])

AC_CONFIG_MACRO_DIRS([m4])

ACLOCAL="$ACLOCAL $ACLOCAL_FLAGS"


# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_MAKE_SET

# Check Compiler
AC_C_INLINE
AC_C_CONST


# Threading 
AC_OPENMP

# Checks for libraries.
AC_CHECK_LIB([fftw3],  [fftw_plan_dft_r2c_3d], [fftw_d=true],[fftw_d=false; AC_MSG_WARN([used to build double precision tools and library])])
AC_CHECK_LIB([fftw3f], [fftwf_plan_dft_r2c_3d],[fftw_f=true],[fftw_f=false; AC_MSG_WARN([used to build single precision tools and library])])

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



# Checks compiler flags
AC_LANG([C])
DEFAULT_CFLAGS=""
AX_CHECK_COMPILE_FLAG([-pedantic],     [OPT_CFLAGS="$DEFAULT_CFLAGS -pdantic"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-Wall],         [OPT_CFLAGS="$DEFAULT_CFLAGS -Wall"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-std=c99],      [OPT_CFLAGS="$DEFAULT_CFLAGS -std=c99"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-O2],           [OPT_CFLAGS="$DEFAULT_CFLAGS -O2"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-O3],           [OPT_CFLAGS="$DEFAULT_CFLAGS -O3"],[],[-Werror])
#AX_CHECK_COMPILE_FLAG([-Os],           [OPT_CFLAGS="$OPT_CFLAGS -finline-functions"],[],[-Werror])
AC_SUBST(DEFAULT_CFLAGS)

# Checks additional debug floags
AC_LANG([C])
DEBUG_CFLAGS=""
AX_CHECK_COMPILE_FLAG([-ggdb],[DEBUG_CFLAGS="$DEBUG_CFLAGS -ggdb"], [],[-Werror])
AX_CHECK_COMPILE_FLAG([-O0],  [DEBUG_CFLAGS="$DEBUG_CFLAGS -O0"],   [],[-Werror])
AC_SUBST(DEBUG_CFLAGS)


# Checks optimization
AC_LANG([C])
OPT_CFLAGS=""
AX_CHECK_COMPILE_FLAG([-finline-functions],           [OPT_CFLAGS="$OPT_CFLAGS -finline-functions"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-funroll-loops],               [OPT_CFLAGS="$OPT_CFLAGS -funroll-loops"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-fprefetch-loop-arrays],       [OPT_CFLAGS="$OPT_CFLAGS -fprefetch-loop-arrays"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-ftree-vectorize],             [OPT_CFLAGS="$OPT_CFLAGS -ftree-vectorize"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-funswitch-loops],             [OPT_CFLAGS="$OPT_CFLAGS -funswitch-loops"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-fgcse-after-reload],          [OPT_CFLAGS="$OPT_CFLAGS -fgcse-after-reload"],[],[-Werror])

AX_CHECK_COMPILE_FLAG([-fno-math-errno],              [OPT_CFLAGS="$OPT_CFLAGS -fno-math-errno"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-fno-trapping-math],           [OPT_CFLAGS="$OPT_CFLAGS -fno-trapping-math"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-ffinite-math-only],           [OPT_CFLAGS="$OPT_CFLAGS -ffinite-math-only"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-fno-rounding-math],           [OPT_CFLAGS="$OPT_CFLAGS -fno-rounding-math"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-fno-signaling-nans],          [OPT_CFLAGS="$OPT_CFLAGS -fno-signaling-nans"],[],[-Werror])

AX_CHECK_COMPILE_FLAG([-funsafe-loop-optimizations],  [OPT_CFLAGS="$OPT_CFLAGS -funsafe-loop-optimizations"],[],[-Werror])
AX_CHECK_COMPILE_FLAG([-Wunsafe-loop-optimizations],  [OPT_CFLAGS="$OPT_CFLAGS -Wunsafe-loop-optimizations"],[],[-Werror])

AX_CHECK_COMPILE_FLAG([-funsafe-math-optimizations ], [OPT_CFLAGS="$OPT_CFLAGS -funsafe-math-optimizations"],[],[-Werror])
AC_SUBST(OPT_CFLAGS)


# Debug
AC_ARG_ENABLE([debug],[  --enable-debug          turn on debugging],
	[case "${enableval}" in
  		yes) debug=true ;;
    		no)  debug=false ;;
      		*) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
      	esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test "x$debug" = xtrue])

# Verbose
AC_ARG_ENABLE([verbose],[  --enable-verbose        turn on verbose mode],
	[case "${enableval}" in
  		yes) verbose=true ;;
    		no)  vebose=false ;;
      		*) AC_MSG_ERROR([bad value ${enableval} for --enable-verbose]) ;;
      	esac],[verbose=false])
AM_CONDITIONAL([VERBOSE], [test "x$verbose" = xtrue])


# Precision
AC_ARG_ENABLE([single],[  --enable-single         turn on single precision],
	[case "${enableval}" in
  		yes) single=true ;;
    		no)  single=false ;;
      		*) AC_MSG_ERROR([bad value ${enableval} for --enable-single]) ;;
      	esac],[single=false])
AM_CONDITIONAL([SINGLE], [test "x$single" = xtrue])


AC_CONFIG_FILES([Makefile src/Makefile])


if [ [ "x$single" = "xtrue" ] && [ "x$fftw_f" = "xfalse" ] ];
then
	AC_MSG_ERROR(["To build single precision tools and lib fftw3 is missing: see README"])
else 
	if [ [ "x$fftw_d" = "xfalse" ] ];then
		AC_MSG_ERROR(["To build double precision tools and lib fftw3 is missing: see README"])
	fi
fi

AC_OUTPUT
