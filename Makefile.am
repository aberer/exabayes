AUTOMAKE_OPTIONS=foreign


WARN_FLAGS= --pedantic -Wall #  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast -Wnested-externs  -Wold-style-definition -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast #   -Wstrict-prototypes   -Wdeclaration-after-statement  -Wmissing-prototypes  
LANG_FLAGS=-std=gnu99 # because of asm 


# TODO they cannot be serious with that...this is the only way of getting rid of the default -O2 -g  
CFLAGS=
CXXFLAGS=

AM_CFLAGS=-D__SIM_SSE3 -D_BAYESIAN -msse3 -D_GNU_SOURCE  -D_OPTIMIZED_FUNCTIONS

if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops
endif

# TODO avx, threads
  # src/pll/avxLikelihood-pll.c src/pll/axml-pll.c 


srcpre=./src
srcxpre=./src/examl
srcppre=./src/pll

SRC=$(srcpre)/nclConfigReader.cpp $(srcpre)/configParser.c $(srcpre)/randomness.c $(srcpre)/bayes.c $(srcpre)/output.c
bin_PROGRAMS=exabayes


if WITH_PLL 
########
# PLL  #
########
OPT_SRC= $(srcpre)/main-pll.c $(srcpre)/main-common.c  $(srcppre)/evaluateGenericSpecial-pll.c $(srcppre)/fastDNAparsimony-pll.c $(srcppre)/makenewzGenericSpecial-pll.c  $(srcppre)/models-pll.c $(srcppre)/optimizeModel-pll.c  $(srcppre)/recom-pll.c $(srcppre)/searchAlgo-pll.c  $(srcppre)/trash-pll.c   $(srcppre)/utils-pll.c  $(srcppre)/bipartitionList-pll.c  $(srcppre)/evaluatePartialGenericSpecial-pll.c  $(srcppre)/genericParallelization-pll.c  $(srcppre)/mem_alloc-pll.c $(srcppre)/newviewGenericSpecial-pll.c  $(srcppre)/randomTree-pll.c $(srcppre)/restartHashTable-pll.c  $(srcppre)/topologies-pll.c  $(srcppre)/treeIO-pll.c 
OPT_HEAD= -I$(srcppre)
FEATURE_FLAGS= -D_USE_PTHREADS -pthread #  -D_USE_ALLREDUCE
else  
#########
# EXAML #
#########
CC=mpicc
CXX=mpicxx
OPT_SRC= $(srcpre)/main-examl.c  $(srcpre)/main-common.c  $(srcxpre)/bipartitionList.c $(srcxpre)/evaluatePartialGenericSpecial.c $(srcxpre)/models.c $(srcxpre)/optimizeModel.c $(srcxpre)/searchAlgo.c $(srcxpre)/trash.c \
 $(srcxpre)/evaluateGenericSpecial.c $(srcxpre)/makenewzGenericSpecial.c $(srcxpre)/newviewGenericSpecial.c $(srcxpre)/restartHashTable.c $(srcxpre)/topologies.c $(srcxpre)/treeIO.c
OPT_HEAD= -I$(srcxpre)
FEATURE_FLAGS= -D_USE_ALLREDUCE
endif 



exabayes_SOURCES= $(SRC) $(OPT_SRC)

AM_CFLAGS+=$(OPT_HEAD)
AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS) $(LANG_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)
AM_LDFLAGS = -lm   -lncl  -L/usr/local/lib/ncl/
