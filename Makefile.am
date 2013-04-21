AUTOMAKE_OPTIONS=foreign
SUBDIRS = ./lib/ncl-2.1.18/

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

bin_PROGRAMS= exabayes #     prsf #  asdsf

include src-helper.am
include src.am 

.PHONY: my-update-src doc 

noinst_LIBRARIES=libbayes.a 



exabayes_LDADD=

# configure examl/pll specific stuff 
if WITH_PLL 
OPT_HEAD= -I./src/pll
libpll_a_SOURCES= ./src/pll/evaluateGenericSpecial-pll.c ./src/pll/axml-pll.c ./src/pll/trash-pll.c ./src/pll/vec_unit.h ./src/pll/makenewzGenericSpecial-pll.c ./src/pll/mem_alloc-pll.c ./src/pll/restartHashTable-pll.c ./src/pll/axml-pll.h ./src/pll/mem_alloc.h ./src/pll/recom-pll.c ./src/pll/searchAlgo-pll.c ./src/pll/phylip_parser/msa_sites.h ./src/pll/phylip_parser/xalloc.h ./src/pll/phylip_parser/phylip.h ./src/pll/phylip_parser/lexer.h ./src/pll/phylip_parser/ssort.h ./src/pll/utils-pll.c ./src/pll/treeIO-pll.c ./src/pll/cycle.h ./src/pll/aligned_buffer.h ./src/pll/globalVariables.h ./src/pll/topologies-pll.c ./src/pll/utils.h ./src/pll/evaluatePartialGenericSpecial-pll.c ./src/pll/fastDNAparsimony-pll.c ./src/pll/randomTree-pll.c ./src/pll/optimizeModel-pll.c ./src/pll/treeIO.h ./src/pll/newviewGenericSpecial-pll.c ./src/pll/bipartitionList-pll.c ./src/pll/models-pll.c ./src/pll/genericParallelization.h 
# FEATURE_FLAGS= -D_USE_PTHREADS
noinst_LIBRARIES+=libpll.a
exabayes_LDADD+=libpll.a
libpll_a_CFLAGS=  -w  $(AM_CFLAGS) $(AM_CXXFLAGS)
else
CC=mpicc -std=c99
CXX=mpicxx 
OPT_HEAD= -Isrc/examl/
FEATURE_FLAGS= -D_USE_ALLREDUCE
noinst_LIBRARIES+= libexaml.a
exabayes_LDADD+=libexaml.a
libexaml_a_CFLAGS= -w  $(AM_CFLAGS) $(AM_CXXFLAGS)   -fgnu89-inline
libexaml_a_SOURCES= ./src/examl/searchAlgo.c ./src/examl/treeIO.c ./src/examl/restartHashTable.c ./src/examl/optimizeModel.c ./src/examl/newviewGenericSpecial.c ./src/examl/trash.c ./src/examl/topologies.c ./src/examl/axml-examl.h ./src/examl/globalVariables.h ./src/examl/evaluatePartialGenericSpecial.c ./src/examl/bipartitionList.c ./src/examl/models.c ./src/examl/makenewzGenericSpecial.c ./src/examl/evaluateGenericSpecial.c ./src/examl/axml.c
endif 

exabayes_LDADD+= libbayes.a ./lib/ncl-2.1.18/ncl/.libs/libncl.a 
exabayes_SOURCES = ./src/exabayes.cpp

AM_CFLAGS+=$(OPT_HEAD) -Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)  -std=c++0x


AM_LDFLAGS=-lm

EXTRA_DIST=lib/Random123-1.07 src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c src/pll/axml-pll.h src/examl/axml-examl.h doc/*



################################################################
# PRSF

prsf_LDADD= libpll.a libbayes.a
prsf_SOURCES=./src/prsf.h ./src/prsf.cpp



################################################################
# ASDSF


# asdsf_SOURCES=src/AvgSplitFreqAssessor.cpp src/AvgSplitFreqAssessor.hpp
# asdsf_CXXFLAGS=$(AM_CXXFLAGS) -D_ASDSF_STANDALONE 
# asdsf_LDADD=./lib/ncl-2.1.18/ncl/.libs/libncl.a  


# if WITH_PLL	
#  asdsf_SOURCES+=./src/pll/mem_alloc-pll.c ./src/pll/bipartitionList-pll.c
# else 
#  asdsf_SOURCES+=./src/examl/bipartitionList.c
# endif





################################################################
# cleanup 

clean:
	-rm -rf *.o *.a  #   $(exabayes_OBJECTS) $(libbayes_A_OBJECTS) libbayes.a $(bin_PROGRAMS )  $(prsf_OBJECTS)  

clean-all : 
	-rm -rf $(exabayes_OBJECTS) $(bin_PROGRAMS )  $(prsf_OBJECTS)  *.a  *.o 


doc : 
	doxygen 
