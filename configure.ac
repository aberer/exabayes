#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_REREQ([2.67])
AC_INIT([exabayes],[1.2],[exabayes-at-googlegroups-dot-com]) 

m4_include([./m4/ax_cxx_compile_stdcxx_11.m4])

AC_CONFIG_SRCDIR([src/])

dnl todo-ish 
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([src/config.h])


AC_DEFINE([PROGRAM_NAME], ["ExaBayes"], ["name used for log files"])

AM_SILENT_RULES([yes])


# Checks for programs.

OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
AC_PROG_CC
AC_PROG_CXX
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS


AX_CXX_COMPILE_STDCXX_11
AC_PROG_CC_STDC
AC_PROG_CC_C99
AC_PROG_RANLIB

dnl NOTICE if you have trouble with that, then install autoconf-archive (or similar)
AX_EXT

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




dnl disallow to use sse3 
AC_ARG_ENABLE( [sse], 
	       [  --disable-sse    compiles ExaBayes without any vectorization (slower), even if available], 
 	       [ case $enableval in 
	       no ) GUESS_SSE="yes" ;; 
 	       yes ) GUESS_SSE="no" ;; 
	       * )  GUESS_SSE="no"
 	       esac ] ,
	       [GUESS_SSE="no"]
 )
AM_CONDITIONAL([NO_SSE], [ test "x$GUESS_SSE" = "xyes" ])

dnl allow to disable avx 
AC_ARG_ENABLE( [avx],
	       [  --disable-avx   compiles ExaBayes with SSE support (slower) only, altough AVX is available ], 
[ case $enableval in 
  no ) GUESSAVX="yes" ;; 
  yes ) GUESSAVX="no" ;;
  * ) GUESSAVX="no"
  esac ] , 
  [GUESSAVX="no"]
  )
AM_CONDITIONAL([NO_AVX], [ test "x$GUESSAVX" = "xyes" ])

AC_ARG_ENABLE( [mpi], 	      
	   [  --enable-mpi  compiles ExaBayes with MPI support (please provide mpi-compiler wrappers to CC/CXX) ],
[  case "${enableval}" in 
   no )  CODE="no" ;; 
   yes ) CODE="yes" ;; 
   *) CODE="no" 
   esac ] ,
[CODE="no"]
)


dnl test, if we have PLL AND force the user to provide a mpi-comiler, if we dont
if test "x$CODE" = "xyes" ; then 
   AC_DEFINE([HAVE_PLL], [0], "")

   AC_CHECK_HEADER(mpi.h,
    [],
    [AC_MSG_ERROR([ 

       Could not find mpi.h. If you build ExaBayes with MPI support, please provide compiler wrappers, e.g., 
              ./configure --enable-mpi CC=mpicc CXX=mpicxx 
       , where specifically mpicxx could also be mpiCC or mpic++ (depends on your installation). 
    ])]
    )

   cache_var=AS_TR_SH([ac_cv_header_mpi.h])
   $as_unset $cache_var

   AC_LANG_PUSH([C++])
   AC_CHECK_HEADER(mpi.h,
    [],
    [AC_MSG_ERROR([ 

       Could not find mpi.h. If you build ExaBayes with MPI support, please provide compiler wrappers, e.g., 
              ./configure --enable-mpi CC=mpicc CXX=mpicxx 
       , where specifically mpicxx could also be mpiCC or mpic++ (depends on your installation). 
    ])]
    )
   AC_LANG_POP([C++]) 

elif test "x$CODE" = "xno" ; then 
   AC_DEFINE([HAVE_PLL], [1], "")
else 
   AC_MSG_ERROR(bad value --enable-mpi) 
fi 

AM_CONDITIONAL([WITH_EXAML], [ test "x$CODE" = "xyes" ])


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
