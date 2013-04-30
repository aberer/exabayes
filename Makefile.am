AUTOMAKE_OPTIONS=foreign
SUBDIRS = ./lib/ncl-2.1.18/

WARN_FLAGS=  -Wall --pedantic  -Wno-unused-function -Wno-write-strings -Wno-vla

# SIMD_FLAGS is potentially problematic: auto-avx vectorization was considered harmfull? 
AM_CFLAGS= $(SIMD_FLAGS) -D_GNU_SOURCE -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3 -D_USE_ALLREDUCE
# TODO 

if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops
endif

# TODO threads
  #  $(srcppre)/genericParallelization-pll.c

bin_PROGRAMS= exabayes prsf asdsf 

include src-helper.am
include src.am 

.PHONY: my-update-src doc 

AM_CFLAGS+= -Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)  -std=c++0x
AM_LDFLAGS=-lm


exabayes_LDADD= 
exabayes_SOURCES = ./src/exabayes.cpp ./src/CommandLine.cpp ./src/CommandLine.hpp

PLLINCLUDE=-I./src/pll
libpll_a_SOURCES= ./src/pll/evaluateGenericSpecial-pll.c ./src/pll/axml-pll.c ./src/pll/trash-pll.c ./src/pll/vec_unit.h ./src/pll/makenewzGenericSpecial-pll.c ./src/pll/mem_alloc-pll.c ./src/pll/restartHashTable-pll.c ./src/pll/axml-pll.h ./src/pll/mem_alloc.h ./src/pll/recom-pll.c ./src/pll/searchAlgo-pll.c ./src/pll/phylip_parser/msa_sites.h ./src/pll/phylip_parser/xalloc.h ./src/pll/phylip_parser/phylip.h ./src/pll/phylip_parser/lexer.h ./src/pll/phylip_parser/ssort.h ./src/pll/utils-pll.c ./src/pll/treeIO-pll.c ./src/pll/cycle.h ./src/pll/aligned_buffer.h ./src/pll/globalVariables.h ./src/pll/topologies-pll.c ./src/pll/utils.h ./src/pll/evaluatePartialGenericSpecial-pll.c ./src/pll/fastDNAparsimony-pll.c ./src/pll/randomTree-pll.c ./src/pll/optimizeModel-pll.c ./src/pll/treeIO.h ./src/pll/newviewGenericSpecial-pll.c ./src/pll/bipartitionList-pll.c ./src/pll/models-pll.c ./src/pll/genericParallelization.h ./src/avxConditional.c
libpll_a_CFLAGS= -w $(AM_CFLAGS) $(PLLINCLUDE)

EXAMLINCLUDE=-I./src/examl/
libexaml_a_SOURCES= ./src/examl/searchAlgo.c ./src/examl/treeIO.c ./src/examl/restartHashTable.c ./src/examl/optimizeModel.c ./src/examl/newviewGenericSpecial.c ./src/examl/trash.c ./src/examl/topologies.c ./src/examl/axml-examl.h ./src/examl/globalVariables.h ./src/examl/evaluatePartialGenericSpecial.c ./src/examl/bipartitionList.c ./src/examl/models.c ./src/examl/makenewzGenericSpecial.c ./src/examl/evaluateGenericSpecial.c ./src/examl/axml.c ./src/avxConditional.c
libexaml_a_CFLAGS= -w $(AM_CFLAGS) -fgnu89-inline  $(EXAMLINCLUDE) $(MPI_CFLAGS)

noinst_LIBRARIES=libbayes.a 

# configure examl/pll specific stuff 
if WITH_EXAML 
LIBTOUSE=libexaml.a
AM_CFLAGS += $(EXAMLINCLUDE) $(MPI_CFLAGS) 
AM_LDFLAGS+=$(MPI_CLDFLAGS)
else
LIBTOUSE=libpll.a
AM_CFLAGS+=$(PLLINCLUDE)
endif 

noinst_LIBRARIES+=$(LIBTOUSE)
exabayes_LDADD+=$(LIBTOUSE)

exabayes_LDADD+=libbayes.a ./lib/ncl-2.1.18/ncl/.libs/libncl.a 
EXTRA_DIST=lib/Random123-1.07 src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c src/pll/axml-pll.h src/examl/axml-examl.h doc/*

################################################################
# PRSF

prsf_LDADD = libbayes.a $(LIBTOUSE) 
prsf_SOURCES= ./src/prsf.cpp



################################################################

# ASDSF
asdsf_LDADD=    libbayes.a $(LIBTOUSE)   ./lib/ncl-2.1.18/ncl/.libs/libncl.a
asdsf_SOURCES= ./src/asdsf.cpp

################################################################
# examl parser 
pre1=./lib/ExaML/parser/
parserExaml_SOURCES=$(pre1)/axml.c c  $(pre1)/parsePartitions.c $(pre1)/axml.h  $(pre1)/globalVariables.h


################################################################
# pll parser 
pre2=./lib/phylogenetic-likelihood-library/parser
parserPll_SOURCES=$(pre2)/axml.c $(pre2)/parsePartitions.c $(pre2)/axml.h $(pre2)/globalVariables.h src/pll/mem_alloc.h ./src/pll/mem_alloc-pll.c
parserPll_CFLAGS=  -w -I$(pre2) -I$(PLLINCLUDE) $(AM_CFLAGS) 

################################################################
# cleanup 

clean:
	-rm -rf *.o *.a $(bin_PROGRAMS ) $(exabayes_OBJECTS) $(libbayes_A_OBJECTS) libbayes.a

clean-all : 
	-rm -rf $(bin_PROGRAMS) *.o *.a 


doc : 
	doxygen 
