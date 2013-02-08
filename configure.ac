#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([exabayes],[0.0.1],[bla@bla.org])
AM_INIT_AUTOMAKE(  exabayes , 0.0.1 )
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([src/config.h])


AM_SILENT_RULES([yes])


# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

AX_MPI()

# Checks for libraries.
dnl AC_CHECK_LIB([ncl], [], [], [ AC_MSG_ERROR([Could not find the ncl library]) ])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h stdint.h stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime getcwd getpagesize gettimeofday memset pow sqrt strcasecmp strchr])



AC_ARG_ENABLE(pll,
  [ --enable-pll   enable compilation with PLL instead of ExaML code base (the latter is default)],
  [case "${enableval}" in
     yes | no ) WITH_PLL="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --enable-pll) ;;
   esac],
  [WITH_PLL="no"]
)
AM_CONDITIONAL([WITH_PLL], [test  "x$WITH_PLL" = "xyes"])


AC_ARG_ENABLE(debug,
  [ --enable-debug    compile with debug information and  less optimization ],
  [case "${enableval}" in
     yes | no ) WITH_DEBUG="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
   esac],
  [WITH_DEBUG="no"]
)
AM_CONDITIONAL([WITH_DEBUG], [test  "x$WITH_DEBUG" = "xyes"])


AC_ARG_ENABLE([debug],[ --enable-debug	Build with debug information ][debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue ])


AC_ARG_ENABLE( [examl], 
	       [ --enable-examl	Use the ExaML code base (ExaML is the default setting)],
	       [examl_usage=true] )
AM_CONDITIONAL([ENABLE_EXAML], [test x$examl_usage = xtrue ])


dnl AC_ARG_ENABLE([mpi],
dnl [ --with-mpi	Use the phylogenetic likelihood library instead of the ExaML code base],
dnl [case "${useval}" in 
dnl       yes) mpi_usage=true ;;
dnl       no)  mpi_usage=false ;;
dnl       *) AC_MSG_ERROR([bad value ${useval} for --with-mpi]) ;; 
dnl esac], [mpi_usage=false])
dnl AM_CONDITIONAL([ENABLE_MPI], [test x$mpi_usage = xtrue ])



# ------------------------------------------------------------------------------
# Setup MPI
# ------------------------------------------------------------------------------

dnl AC_ARG_ENABLE(mpi,AS_HELP_STRING([--enable-mpi],[enable mpi code]))
dnl if test "x$enable_mpi" = "xyes" ; then
dnl     AX_MPI((
dnl             AC_DEFINE(_FINE_GRAIN_MPI, [], [Description])
dnl             CFLAGS="$CFLAGS -std=c11"
dnl  AC_MSG_NOTICE([MPI detected])),
dnl            AC_MSG_ERROR([No MPI detected]))
dnl fi
dnl AM_CONDITIONAL([USE_MPI],test "x$enable_mpi" = "xyes")



AC_CONFIG_FILES([Makefile])
AC_OUTPUT
