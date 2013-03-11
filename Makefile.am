AUTOMAKE_OPTIONS=foreign





SUBDIRS = ./lib/ncl-2.1.18/



WARN_FLAGS= --pedantic -Wall #  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast -Wnested-externs  -Wold-style-definition -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast #   -Wstrict-prototypes   -Wdeclaration-after-statement  -Wmissing-prototypes  
LANG_FLAGS=-std=gnu99 # because of asm 


AM_CFLAGS=-D__SIM_SSE3  -msse3 -D_GNU_SOURCE  -D_OPTIMIZED_FUNCTIONS   -D_BAYESIAN

if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops
endif

# TODO avx, threads
  # src/pll/avxLikelihood-pll.c src/pll/axml-pll.c 


# prefixes 
srcpre=./src
srcxpre=./src/examl
srcppre=./src/pll


SRC=$(srcpre)/nclConfigReader.cpp $(srcpre)/configParser.c		\
$(srcpre)/randomness.c $(srcpre)/bayes.c $(srcpre)/output.c		\
$(srcpre)/proposals.c $(srcpre)/chain.c $(srcpre)/convergence.c		\
$(srcpre)/treeRead.c $(srcpre)/bayes-topo.c $(srcpre)/utils-topo.c	\
$(srcpre)/main-common.c 

bin_PROGRAMS=exabayes

# header files
SRC+=src/bayes-topo.h src/common.h src/convergence.h			\
src/main-common.h src/output.h src/proposalStructs.h src/treeRead.h	\
src/chain.h src/globals.h src/nclConfigReader.h src/proposals.h		\
src/randomness.h src/utils-topo.h $(srcpre)/rng.h



if WITH_PLL 
########
# PLL  #
########
OPT_SRC= $(srcpre)/main-pll.c $(srcppre)/evaluateGenericSpecial-pll.c	\
$(srcppre)/fastDNAparsimony-pll.c					\
$(srcppre)/makenewzGenericSpecial-pll.c $(srcppre)/models-pll.c		\
$(srcppre)/optimizeModel-pll.c $(srcppre)/recom-pll.c			\
$(srcppre)/searchAlgo-pll.c $(srcppre)/utils-pll.c			\
$(srcppre)/bipartitionList-pll.c					\
$(srcppre)/evaluatePartialGenericSpecial-pll.c				\
$(srcppre)/mem_alloc-pll.c $(srcppre)/newviewGenericSpecial-pll.c	\
$(srcppre)/randomTree-pll.c $(srcppre)/restartHashTable-pll.c		\
$(srcppre)/treeIO-pll.c $(srcppre)/trash-pll.c				\
$(srcppre)/topologies-pll.c # $(srcppre)/genericParallelization-pll.c

OPT_SRC+= $(srcppre)/aligned_buffer.h $(srcppre)/cycle.h		\
$(srcppre)/globalVariables.h $(srcppre)/treeIO.h			\
$(srcppre)/vec_unit.h $(srcppre)/axml.h					\
$(srcppre)/genericParallelization.h $(srcppre)/mem_alloc.h		\
$(srcppre)/utils.h $(srcppre)/phylip_parser/lexer.h			\
$(srcppre)/phylip_parser/phylip.h $(srcppre)/phylip_parser/xalloc.h	\
$(srcppre)/phylip_parser/msa_sites.h $(srcppre)/phylip_parser/ssort.h

OPT_HEAD= -I$(srcppre)
# FEATURE_FLAGS= -D_USE_PTHREADS -pthread #  -D_USE_ALLREDUCE
else  
#########
# EXAML #
#########
CC=mpicc
CXX=mpicxx
OPT_SRC= $(srcpre)/main-examl.c $(srcpre)/make-random-tree-examl.c	\
 $(srcxpre)/bipartitionList.c						\
 $(srcxpre)/evaluatePartialGenericSpecial.c $(srcxpre)/models.c		\
 $(srcxpre)/optimizeModel.c $(srcxpre)/treeIO.c				\
 $(srcxpre)/searchAlgo.c $(srcxpre)/evaluateGenericSpecial.c		\
 $(srcxpre)/makenewzGenericSpecial.c					\
 $(srcxpre)/newviewGenericSpecial.c $(srcxpre)/restartHashTable.c	\
 $(srcxpre)/trash.c $(srcxpre)/topologies.c
OPT_SRC+=$(srcxpre)/axml.h  src/examl/globalVariables.h
OPT_HEAD= -I$(srcxpre)

FEATURE_FLAGS= -D_USE_ALLREDUCE
endif 



exabayes_SOURCES= $(SRC) $(OPT_SRC)

AM_CFLAGS+=$(OPT_HEAD) -Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS) $(LANG_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)
exabayes_LDADD=./lib/ncl-2.1.18/ncl/.libs/libncl.a
AM_LDFLAGS=-lm


EXTRA_DIST=lib/Random123-1.07 #  ./utils/defaultRun.sh



clean-local: 
	-rm -rf *.o exabayes
