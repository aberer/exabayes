AUTOMAKE_OPTIONS=foreign


WARN_FLAGS=--pedantic -Wall #-Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes  -Wdeclaration-after-statement -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast
LANG_FLAGS=-std=gnu99 # because of asm 


# TODO they cannot be serious with that...this is the only way of getting rid of the default -O2 -g  
CFLAGS=

AM_CFLAGS=-D__SIM_SSE3 -D_BAYESIAN -msse3 -D_GNU_SOURCE  -D_OPTIMIZED_FUNCTIONS

if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops
endif



# TODO avx, threads
  # src/pll/avxLikelihood-pll.c src/pll/axml-pll.c 

# this is fubar
# OPT_SRC=lib/ncl/ncl/nxsassumptionsblock.cpp   lib/ncl/ncl/nxsdatablock.cpp       lib/ncl/ncl/nxspublicblocks.cpp  lib/ncl/ncl/nxstaxaassociationblock.cpp  lib/ncl/ncl/nxsunalignedblock.cpp lib/ncl/ncl/nxsblock.cpp              lib/ncl/ncl/nxsdistancesblock.cpp  lib/ncl/ncl/nxsreader.cpp        lib/ncl/ncl/nxstaxablock.cpp lib/ncl/ncl/nxscharactersblock.cpp    lib/ncl/ncl/nxsexception.cpp       lib/ncl/ncl/nxssetreader.cpp     lib/ncl/ncl/nxstoken.cpp lib/ncl/ncl/nxscxxdiscretematrix.cpp  lib/ncl/ncl/nxsmultiformat.cpp     lib/ncl/ncl/nxsstring.cpp        lib/ncl/ncl/nxstreesblock.cpp
# OPT_HEAD= -Ilib/ncl-2.1.18/ncl



if WITH_PLL 
########
# PLL  #
########
OPT_SRC= src/main-pll.c  src/nclConfigReader.cpp src/configParser.c src/randomness.c src/pll/evaluateGenericSpecial-pll.c src/pll/fastDNAparsimony-pll.c src/pll/makenewzGenericSpecial-pll.c  src/pll/models-pll.c src/pll/optimizeModel-pll.c  src/pll/recom-pll.c src/pll/searchAlgo-pll.c  src/pll/trash-pll.c   src/pll/utils-pll.c  src/pll/bipartitionList-pll.c  src/pll/evaluatePartialGenericSpecial-pll.c  src/pll/genericParallelization-pll.c  src/pll/mem_alloc-pll.c src/pll/newviewGenericSpecial-pll.c  src/pll/randomTree-pll.c src/pll/restartHashTable-pll.c  src/pll/topologies-pll.c  src/pll/treeIO-pll.c 
OPT_HEAD= -I./src/pll/
bin_PROGRAMS=exabayes_pll
exabayes_pll_SOURCES= src/bayes.c $(OPT_SRC)
FEATURE_FLAGS= -D_USE_PTHREADS -pthread #  -D_USE_ALLREDUCE
else  
#########
# EXAML #
#########
CC=mpicc
CXX=mpicxx
OPT_SRC= src/main-examl.c src/nclConfigReader.cpp src/main-common.c src/randomness.c src/configParser.c src/examl/bipartitionList.c src/examl/evaluatePartialGenericSpecial.c src/examl/models.c src/examl/optimizeModel.c src/examl/searchAlgo.c src/examl/trash.c \
src/examl-legacy.c src/examl/evaluateGenericSpecial.c src/examl/makenewzGenericSpecial.c src/examl/newviewGenericSpecial.c src/examl/restartHashTable.c src/examl/topologies.c src/examl/treeIO.c
OPT_HEAD= -I./src/examl/
bin_PROGRAMS=exabayes
exabayes_SOURCES=src/bayes.c $(OPT_SRC)
FEATURE_FLAGS= -D_USE_ALLREDUCE
endif 





AM_CFLAGS+=$(OPT_HEAD)
AM_CFLAGS+= $(FEATURE_FLAGS) $(WARN_FLAGS) $(LANG_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)
AM_LDFLAGS = -lm   -lncl  -L/usr/local/lib/ncl/


