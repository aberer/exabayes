AUTOMAKE_OPTIONS=foreign


SUBDIRS = ./lib/ncl-2.1.18/

WARN_FLAGS= --pedantic -Wall #  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast -Wnested-externs  -Wold-style-definition -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast #   -Wstrict-prototypes   -Wdeclaration-after-statement  -Wmissing-prototypes  
LANG_FLAGS=-std=gnu99 # because of asm 


# SIMD_FLAGS is potentially problematic: auto-avx vectorization was considered harmfull? 
AM_CFLAGS= $(SIMD_FLAGS) -D_GNU_SOURCE -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3 -D_BAYESIAN

if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops
endif

# TODO threads
  #  $(srcppre)/genericParallelization-pll.c


# prefixes 
srcxpre=./src/examl
srcppre=./src/pll


SRC=src/nclConfigReader.cpp src/configParser.c		\
src/randomness.c src/bayes.c src/output.c		\
src/proposals.c src/chain.c src/convergence.c		\
src/treeRead.c src/bayes-topo.c src/utils-topo.c	\
src/main-common.c  src/eval.c src/adapterCode.c src/main.c  src/avxConditional.c

bin_PROGRAMS=exabayes

# header files
SRC+=src/bayes-topo.h src/common.h src/convergence.h			\
 src/main-common.h src/output.h src/proposalStructs.h src/treeRead.h	\
 src/chain.h src/globals.h src/nclConfigReader.h src/proposals.h	\
 src/randomness.h src/utils-topo.h src/rng.h src/eval.h	\
 src/adapterCode.h src/axml.h 



if WITH_PLL 
########
# PLL  #
########
OPT_SRC= $(srcppre)/evaluateGenericSpecial-pll.c	\
$(srcppre)/fastDNAparsimony-pll.c					\
$(srcppre)/makenewzGenericSpecial-pll.c $(srcppre)/models-pll.c		\
$(srcppre)/optimizeModel-pll.c $(srcppre)/recom-pll.c			\
	$(srcppre)/searchAlgo-pll.c			\
$(srcppre)/bipartitionList-pll.c					\
$(srcppre)/evaluatePartialGenericSpecial-pll.c				\
$(srcppre)/mem_alloc-pll.c $(srcppre)/newviewGenericSpecial-pll.c	\
$(srcppre)/randomTree-pll.c $(srcppre)/restartHashTable-pll.c		\
$(srcppre)/treeIO-pll.c $(srcppre)/trash-pll.c				\
$(srcppre)/topologies-pll.c    $(srcppre)/utils-pll.c

OPT_SRC+= $(srcppre)/aligned_buffer.h $(srcppre)/cycle.h		\
$(srcppre)/globalVariables.h $(srcppre)/treeIO.h			\
$(srcppre)/vec_unit.h 					\
$(srcppre)/genericParallelization.h $(srcppre)/mem_alloc.h		\
$(srcppre)/utils.h $(srcppre)/phylip_parser/lexer.h			\
$(srcppre)/phylip_parser/phylip.h $(srcppre)/phylip_parser/xalloc.h	\
$(srcppre)/phylip_parser/msa_sites.h $(srcppre)/phylip_parser/ssort.h

OPT_HEAD= -I$(srcppre)
# FEATURE_FLAGS= -D_USE_PTHREADS
else  
#########
# EXAML #
#########
CC=mpicc
CXX=mpicxx
OPT_SRC=  src/make-random-tree-examl.c	\
 $(srcxpre)/bipartitionList.c						\
 $(srcxpre)/evaluatePartialGenericSpecial.c $(srcxpre)/models.c		\
 $(srcxpre)/optimizeModel.c $(srcxpre)/treeIO.c				\
 $(srcxpre)/searchAlgo.c $(srcxpre)/evaluateGenericSpecial.c		\
 $(srcxpre)/makenewzGenericSpecial.c					\
 $(srcxpre)/newviewGenericSpecial.c $(srcxpre)/restartHashTable.c	\
 $(srcxpre)/trash.c $(srcxpre)/topologies.c
OPT_SRC+=  src/examl/globalVariables.h
OPT_HEAD= -I$(srcxpre)
FEATURE_FLAGS= -D_USE_ALLREDUCE
endif 

exabayes_SOURCES= $(SRC) $(OPT_SRC)

AM_CFLAGS+=$(OPT_HEAD) -Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS) $(LANG_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)
exabayes_LDADD=./lib/ncl-2.1.18/ncl/.libs/libncl.a
AM_LDFLAGS=-lm

EXTRA_DIST=lib/Random123-1.07 src/main-examl.c src/main-pll.c		\
              src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c	\
              src/pll/axml-pll.h src/examl/axml-examl.h


clean-local: 
	-rm -rf *.o exabayes
