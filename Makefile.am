AUTOMAKE_OPTIONS=foreign
SUBDIRS = ./lib/ncl-2.1.18/


# -f may also help with examl  
WARN_FLAGS=  -Wall --pedantic  -Wno-unused-function -Wno-write-strings -Wno-vla

# SIMD_FLAGS is potentially problematic: auto-avx vectorization was considered harmfull? 
AM_CFLAGS= $(SIMD_FLAGS) -D_GNU_SOURCE -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3


if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops 
endif

# TODO threads
  #  $(srcppre)/genericParallelization-pll.c

bin_PROGRAMS= exabayes   prsf 

include src-helper.am
include src.am 

.PHONY: my-update-src doc 


exabayes_LDADD=./lib/ncl-2.1.18/ncl/.libs/libncl.a

# configure examl/pll specific stuff 
if WITH_PLL 
OPT_HEAD= -I./src/pll
# FEATURE_FLAGS= -D_USE_PTHREADS
exabayes_SOURCES+=$(libpll_a_SOURCES) 
else
CC=mpicc -std=c99
CXX=mpicxx 
OPT_HEAD= -Isrc/examl/
FEATURE_FLAGS= -D_USE_ALLREDUCE
exabayes_SOURCES+=$(libexaml_a_SOURCES)
endif 

AM_CFLAGS+=$(OPT_HEAD) -Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)  -std=c++0x


AM_LDFLAGS=-lm

EXTRA_DIST=lib/Random123-1.07 src/main-examl_cond.c src/main-pll_cond.c		\
              src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c	\
              src/pll/axml-pll.h src/examl/axml-examl.h doc/*


if WITH_PLL
prsf_SOURCES+=./src/pll/mem_alloc-pll.c ./src/pllmem_alloc-pll.h
endif

prsf_CXXFLAGS=$(AM_CXXFLAGS) -D_STANDALONE


clean:
	-rm -rf $(exabayes_OBJECTS) $(bin_PROGRAMS )  $(prsf_OBJECTS) 


doc : 
	doxygen 
