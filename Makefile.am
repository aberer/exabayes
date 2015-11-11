AUTOMAKE_OPTIONS= foreign subdir-objects

EXTRA_DIST=./lib/Random123-1.08 ./examples/* ./manual/manual.pdf ./manual/css/org.css ./manual/img/banner.png ./manual/img/banner2.png ./manual/manual.html ./README.txt ./LICENSE.txt ./COPYING ./bin/ ./lib/pll/avxLikelihood.c

noinst_LIBRARIES=libbayes.a libncl.a libcommdummy.a

.PHONY: my-update-src doc mydist  man  dist

WARN_FLAGS= -Wno-unused-variable -Wall --pedantic -Wextra  -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable #    -Wuninitialized  -Wno-mismatched-tags -Wno-unused-private-field
RAX_FEATURE=-D_GNU_SOURCE -fomit-frame-pointer -funroll-loops -D__STRICT_ANSI__ 

AM_CXXFLAGS=$(RAX_FEATURE) $(WARN_FLAGS) $(PTHREAD_CFLAGS)
AM_CFLAGS=$(RAX_FEATURE) $(WARN_FLAGS) -std=c1x
AM_LDFLAGS=$(PTHREAD_LDFLAGS) $(AM_LDLFLAGS)

if NO_SSE
VECTOR_FLAGS=-DMANUAL_SSE_OVERRIDE -DMANUAL_AVX_OVERRIDE -mno-avx -mno-sse3 # $(SIMD_FLAGS) -march=native -mtune=generic -mno-avx -mno-sse3
else 
if NO_AVX
VECTOR_FLAGS=-DMANUAL_AVX_OVERRIDE -D_OPTIMIZED_FUNCTIONS -msse3 -mno-avx  # $(SIMD_FLAGS) -march=native -mtune=native -mno-avx -m
else 
VECTOR_FLAGS= -D_OPTIMIZED_FUNCTIONS $(SIMD_FLAGS) -march=native -mtune=native
endif
endif

SRCINCLUDE=-I$(top_srcdir)/src

BAYES_INCLUDES=-I$(top_srcdir)/lib/Random123-1.08/include -I$(top_srcdir)/lib/ncl-2.1.18/ $(SRCINCLUDE) -I$(top_srcdir)/src/comm/dummy
BAYES_INCLUDES_PARA=-I$(top_srcdir)/lib/Random123-1.08/include -I$(top_srcdir)/lib/ncl-2.1.18/  $(SRCINCLUDE) -I$(top_srcdir)/src/comm/mpi

include src-helper.am
include src.am 

# the comm stuff 
libcommdummy_a_CXXFLAGS=-O3 $(AM_CXXFLAGS) $(VECTOR_FLAGS) $(BAYES_INCLUDES) 

# the new pll 
libpll_a_CFLAGS=$(AM_CFLAGS) -w $(VECTOR_FLAGS) -O2 $(SRCINCLUDE)  
noinst_LIBRARIES+=libpll.a

bin_PROGRAMS = yggdrasil  sdsf postProcParam credibleSet extractBips consense parser 	

# configure examl/pll specific stuff 
if WITH_MPI 
bin_PROGRAMS+=exabayes 
noinst_LIBRARIES+=libcommmpi.a
libcommmpi_a_CXXFLAGS=$(AM_CXXFLAGS) -O3 $(VECTOR_FLAGS) $(MPI_CXXFLAGS) $(BAYES_INCLUDES_PARA) 
endif 

if WITH_TEST
 bin_PROGRAMS+=exatest
endif 

seqlibs= libncl.a libbayes.a libcommdummy.a  libpll.a libncl.a  
paralibs= libbayes.a libcommmpi.a libpll.a libncl.a

libbayes_a_CXXFLAGS= $(AM_CXXFLAGS) -O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS)    
libbayes_a_CFLAGS=-O2 $(BAYES_INCLUDES) $(VECTOR_FLAGS)  

################################################################
# exabayes 
exabayes_LDADD= $(paralibs) $(MPI_CXXLDFLAGS)
exabayes_SOURCES = ./src/apps/exabayes.cpp
exabayes_CXXFLAGS= $(AM_CXXFLAGS) -O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS)   
exabayes_LDFLAGS=$(MPI_CXXLDFLAGS)

################################################################
# yggdrasil 
yggdrasil_LDADD=$(seqlibs)
yggdrasil_SOURCES = ./src/apps/exabayes.cpp
yggdrasil_CXXFLAGS=$(AM_CXXFLAGS) -O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS) -D_IS_YGG

################################################################
# SDSF 
sdsf_CXXFLAGS= $(AM_CXXFLAGS) -O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS) 
sdsf_LDADD= libncl.a libbayes.a    libpll.a libncl.a 
sdsf_SOURCES= ./src/apps/sdsf.cpp 

################################################################
# the parser 
parser_SOURCES=./src/apps/parser.cpp
parser_LDADD=$(seqlibs)
parser_CXXFLAGS=-O3 $(AM_CXXFLAGS) $(BAYES_INCLUDES) $(VECTOR_FLAGS)  -Wno-deprecated-writable-strings  

################################################################
# extractBips 
extractBips_CXXFLAGS=-O3 $(AM_CXXFLAGS) $(VECTOR_FLAGS) $(BAYES_INCLUDES)  
extractBips_LDADD=$(seqlibs)
extractBips_SOURCES=./src/apps/extractBips.cpp  

################################################################
# consense 
consense_CXXFLAGS=-O3 $(AM_CXXFLAGS) $(VECTOR_FLAGS) $(BAYES_INCLUDES) 
consense_LDADD=$(seqlibs)
consense_SOURCES=./src/apps/consense.cpp

################################################################
# credibleSet 
credibleSet_CXXFLAGS=-O3 $(AM_CXXFLAGS) $(VECTOR_FLAGS) $(BAYES_INCLUDES) 
credibleSet_LDADD= $(seqlibs)
credibleSet_SOURCES=./src/apps/credibleSet.cpp

################################################################
# postProcParam 
postProcParam_CXXFLAGS=-O3 $(AM_CXXFLAGS) $(VECTOR_FLAGS) $(BAYES_INCLUDES) 
postProcParam_LDADD=$(seqlibs)
postProcParam_SOURCES=./src/apps/postProcParam.cpp

################################################################
# ncl 
libncl_a_CXXFLAGS=-w $(AM_CXXFLAGS) -O3 -I$(top_srcdir)/lib/ncl-2.1.18/ -D__STRICT_ANSI__   

################################################################

# test 
exatest_SOURCES = ./tests/exatest.cpp
exatest_CXXFLAGS=$(AM_CXXFLAGS) -O3 $(VECTOR_FLAGS) $(BAYES_INCLUDES)  -D_IS_YGG
exatest_LDFLAGS=-L/usr/local/lib -lgtest -lpthread 
exatest_LDADD= $(seqlibs) 

if WITH_MPI
 exatest_CXXFLAGS+=-D_WITH_MPI	
 # teststr="mpirun -np 2 ./exatest"
 teststr="./exatest"
else 
 teststr="./exatest"
endif

test : exatest
	$(RM) ExaBayes_*test-*
	eval $(teststr)

################################################################

doc : 
	doxygen 

man:  
	$(top_srcdir)/manual/process.sh $(top_srcdir)/manual

date :
	$(top_srcdir)/utils/updateDate.sh

dist-hook: 
	make date 
