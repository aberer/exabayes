AUTOMAKE_OPTIONS=foreign
SUBDIRS = ./lib/ncl-2.1.18/


# -f may also help with examl  
WARN_FLAGS=  -Wall --pedantic -fgnu89-inline  -Wno-unused-function

# SIMD_FLAGS is potentially problematic: auto-avx vectorization was considered harmfull? 
AM_CFLAGS= $(SIMD_FLAGS) -D_GNU_SOURCE -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3 -D_BAYESIAN


if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops
endif

# TODO threads
  #  $(srcppre)/genericParallelization-pll.c

bin_PROGRAMS=exabayes

include src-helper.am
include src.am 

.PHONY: my-update-src doc 


exabayes_LDADD=./lib/ncl-2.1.18/ncl/.libs/libncl.a

# configure examl/pll specific stuff 
if WITH_PLL 
OPT_HEAD= -I./src/pll
# FEATURE_FLAGS= -D_USE_PTHREADS
exabayes_LDADD+=libpll.a
noinst_LIBRARIES=libpll.a
# exabayes_SOURCES+= src/make-random-tree_cond.c
else
CC=mpicc -std=gnu99
CXX=mpic++
OPT_HEAD= -Isrc/examl/
FEATURE_FLAGS= -D_USE_ALLREDUCE
# exabayes_SOURCES+=$(EXAML_SRC) 
exabayes_LDADD+=libexaml.a
noinst_LIBRARIES=libexaml.a
endif 

AM_CFLAGS+=$(OPT_HEAD) -Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS) $(LANG_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)

libpll_a_CFLAGS= -w $(AM_CFLAGS)
libexaml_a_CFLAGS=-w $(AM_CFLAGS)


AM_LDFLAGS=-lm

EXTRA_DIST=lib/Random123-1.07 src/main-examl_cond.c src/main-pll_cond.c		\
              src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c	\
              src/pll/axml-pll.h src/examl/axml-examl.h


clean: 
	-rm -rf *.o exabayes *.a 

clean-exa: 
	-rm -rf $(exabayes_OBJECTS) exabayes 


doc : 
	doxygen 
