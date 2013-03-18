#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_REREQ([2.67])
AC_INIT([exabayes],[0.1],[exabayes@googlegroups.com]) # correct? 
AM_INIT_AUTOMAKE(  exabayes , 0.1 )
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([src/config.h])

AC_DEFINE([PROGRAM_NAME], ["ExaBayes"], ["name used for log files"])

AM_SILENT_RULES([yes])

# configures the ncl config reader 
AC_CONFIG_SUBDIRS([lib/ncl-2.1.18])

# Checks for programs.

OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
AC_PROG_CC
AC_PROG_CXX
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS

AC_PROG_CC_STDC
AC_PROG_CC_C99
AC_PROG_RANLIB


dnl NOTICE if you have trouble with that, then install autoconf-archive (or similar)
AX_EXT

# Checks for libraries.

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


AC_ARG_ENABLE(pll, [--enable-pll   enable compilation with PLL instead of ExaML code base] )
AC_DEFINE([HAVE_PLL], [0], ["is 1, if built with the PLL "])
AS_IF([ test x$enable_pll = xyes ], [AC_DEFINE(HAVE_PLL,1)]) 
AM_CONDITIONAL([WITH_PLL], [ test x$enable_pll = xyes ] )

AC_ARG_ENABLE(mydebug,[ --enable-mydebug	Build with debug information ])
AM_CONDITIONAL([WITH_DEBUG], [ test x$enable_mydebug = xyes ])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
