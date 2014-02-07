AUTOMAKE_OPTIONS= foreign subdir-objects

EXTRA_DIST=./lib/Random123-1.08 ./src/pll/avxLikelihood-pll.c ./src/examl/avxLikelihood.c ./src/pll/axml-pll.h ./src/examl/axml-examl.h ./examples/* ./build.sh ./apple-build.sh ./manual/manual.pdf ./manual/css/org.css ./manual/img/banner.png ./manual/manual.html ./README.txt ./LICENSE.txt ./COPYING ./bin/

noinst_LIBRARIES=libbayes.a libncl.a 

.PHONY: my-update-src doc mydist  man  dist

WARN_FLAGS= -Wall --pedantic -Wextra -Wno-unused-function -Wno-unused-parameter  -Wuninitialized -Wno-mismatched-tags -Wno-unused-private-field -Wsign-compare
RAX_FEATURE= -D_GNU_SOURCE -fomit-frame-pointer -funroll-loops -D__STRICT_ANSI__ 

if NO_SSE
VECTOR_FLAGS=-DMANUAL_SSE_OVERRIDE -DMANUAL_AVX_OVERRIDE -mno-avx -mno-sse3   -march=native -mtune=generic 
else 
if NO_AVX
VECTOR_FLAGS=-DMANUAL_AVX_OVERRIDE -mno-avx -D_OPTIMIZED_FUNCTIONS -march=native -mtune=native
else 
VECTOR_FLAGS= -D_OPTIMIZED_FUNCTIONS  -march=native -mtune=native $(SIMD_FLAGS)
endif
endif

RAXML_OPT=-O2

SRCINCLUDE=-I$(top_srcdir)/src

BAYES_INCLUDES=-I$(top_srcdir)/lib/Random123-1.08/include -I$(top_srcdir)/lib/ncl-2.1.18/  $(SRCINCLUDE) 

include src-helper.am
include src.am 

libpll_a_SOURCES= ./src/pll/evaluateGenericSpecial-pll.c ./src/pll/axml-pll.c ./src/pll/trash-pll.c ./src/pll/vec_unit.h ./src/pll/makenewzGenericSpecial-pll.c ./src/pll/mem_alloc-pll.c ./src/pll/restartHashTable-pll.c ./src/pll/axml-pll.h ./src/pll/mem_alloc.h ./src/pll/recom-pll.c ./src/pll/searchAlgo-pll.c ./src/pll/phylip_parser/msa_sites.h ./src/pll/phylip_parser/xalloc.h ./src/pll/phylip_parser/phylip.h ./src/pll/phylip_parser/lexer.h ./src/pll/phylip_parser/ssort.h ./src/pll/utils-pll.c ./src/pll/treeIO-pll.c ./src/pll/cycle.h ./src/pll/aligned_buffer.h ./src/pll/globalVariables.h ./src/pll/topologies-pll.c ./src/pll/utils.h ./src/pll/evaluatePartialGenericSpecial-pll.c ./src/pll/fastDNAparsimony-pll.c ./src/pll/randomTree-pll.c ./src/pll/optimizeModel-pll.c ./src/pll/treeIO.h ./src/pll/newviewGenericSpecial-pll.c ./src/pll/bipartitionList-pll.c ./src/pll/models-pll.c ./src/pll/genericParallelization.h src/avxConditional.c
libpll_a_CFLAGS= -w $(VECTOR_FLAGS) $(RAX_FEATURE) $(RAXML_OPT) $(SRCINCLUDE)  $(AM_CPP_FLAGS)
libexaml_a_SOURCES= ./src/examl/searchAlgo.c ./src/examl/treeIO.c ./src/examl/restartHashTable.c ./src/examl/optimizeModel.c ./src/examl/newviewGenericSpecial.c ./src/examl/trash.c ./src/examl/topologies.c ./src/examl/axml-examl.h ./src/examl/globalVariables.h ./src/examl/evaluatePartialGenericSpecial.c ./src/examl/bipartitionList.c ./src/examl/models.c ./src/examl/makenewzGenericSpecial.c ./src/examl/evaluateGenericSpecial.c ./src/examl/axml.c ./src/avxConditional.c ./src/examl/fastDNAparsimony.c ./src/examl/mpiMock.c ./src/examl/mpiMock.h
libexaml_a_CFLAGS= -w $(VECTOR_FLAGS)  -fgnu89-inline   $(RAX_FEATURE) $(RAXML_OPT)  $(SRCINCLUDE)  $(AM_CPP_FLAGS)

# configure examl/pll specific stuff 
if WITH_EXAML 
LIBTOUSE=libexaml.a
noinst_LIBRARIES+=libexaml.a
bin_PROGRAMS = exabayes 
else
LIBTOUSE=libpll.a
bin_PROGRAMS = yggdrasil postProcParam credibleSet extractBips consense parser asdsf
noinst_LIBRARIES+=libpll.a 
endif 

libbayes_a_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)
libbayes_a_CFLAGS=$(RAXML_OPT)  $(BAYES_INCLUDES) $(VECTOR_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)
exabayes_LDADD= libncl.a libbayes.a 
exabayes_LDADD+=$(LIBTOUSE)
exabayes_SOURCES = ./src/apps/exabayes.cpp ./src/config/CommandLine.cpp ./src/config/CommandLine.hpp
exabayes_CXXFLAGS=-O3  $(BAYES_INCLUDES) $(VECTOR_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)


yggdrasil_LDADD= libncl.a libbayes.a 
yggdrasil_LDADD+=$(LIBTOUSE)
yggdrasil_SOURCES = ./src/apps/exabayes.cpp ./src/config/CommandLine.cpp ./src/config/CommandLine.hpp
yggdrasil_CXXFLAGS=-O3  $(BAYES_INCLUDES) $(VECTOR_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)


if WITH_EXAML
cleanobjs=$(exabayes_OBJECTS)
else 
cleanobjs=$(yggdrasil_OBJECTS)
endif


################################################################
# ASDSF 
asdsf_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS) $(RAX_FEATURE) $(AM_CPP_FLAGS)
asdsf_LDADD= libbayes.a  $(LIBTOUSE)  libncl.a
asdsf_SOURCES= ./src/apps/asdsf.cpp 

################################################################
# the parser 
parser_SOURCES=./src/apps/parser.cpp
parser_LDADD= libbayes.a  $(LIBTOUSE) libncl.a
parser_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(VECTOR_FLAGS) $(RAX_FEATURE) -Wno-deprecated-writable-strings  $(AM_CPP_FLAGS)

################################################################
# extractBips 
extractBips_CXXFLAGS=-O3 $(VECTOR_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
extractBips_LDADD= libbayes.a $(LIBTOUSE)  libncl.a
extractBips_SOURCES=./src/apps/extractBips.cpp  


################################################################
# consense 
consense_CXXFLAGS=-O3 $(VECTOR_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
consense_LDADD= libbayes.a  $(LIBTOUSE)  libncl.a
consense_SOURCES=./src/apps/consense.cpp


################################################################
# credibleSet 
credibleSet_CXXFLAGS=-O3 $(VECTOR_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
credibleSet_LDADD= libbayes.a  $(LIBTOUSE) libncl.a
credibleSet_SOURCES=./src/apps/credibleSet.cpp

################################################################
# postProcParam 
postProcParam_CXXFLAGS=-O3 $(VECTOR_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
postProcParam_LDADD= libbayes.a  $(LIBTOUSE) 
postProcParam_SOURCES=./src/apps/postProcParam.cpp

################################################################
# ncl 
libncl_a_SOURCES = ./lib/ncl-2.1.18/ncl/nxsreader.cpp ./lib/ncl-2.1.18/ncl/nxscxxdiscretematrix.cpp ./lib/ncl-2.1.18/ncl/nxstaxaassociationblock.cpp ./lib/ncl-2.1.18/ncl/nxsunalignedblock.cpp ./lib/ncl-2.1.18/ncl/nxsmultiformat.cpp ./lib/ncl-2.1.18/ncl/nxsdatablock.cpp ./lib/ncl-2.1.18/ncl/nxspublicblocks.cpp ./lib/ncl-2.1.18/ncl/nxstaxablock.cpp ./lib/ncl-2.1.18/ncl/nxscharactersblock.cpp ./lib/ncl-2.1.18/ncl/nxsstring.cpp ./lib/ncl-2.1.18/ncl/nxstoken.cpp ./lib/ncl-2.1.18/ncl/nxsassumptionsblock.cpp ./lib/ncl-2.1.18/ncl/nxsblock.cpp ./lib/ncl-2.1.18/ncl/nxsdistancesblock.cpp ./lib/ncl-2.1.18/ncl/nxssetreader.cpp ./lib/ncl-2.1.18/ncl/nxsexception.cpp ./lib/ncl-2.1.18/ncl/nxstreesblock.cpp ./lib/ncl-2.1.18/ncl/ncl.h ./lib/ncl-2.1.18/ncl/nxscharactersblock.h ./lib/ncl-2.1.18/ncl/nxsdistancedatum.h ./lib/ncl-2.1.18/ncl/nxsreader.h ./lib/ncl-2.1.18/ncl/nxstoken.h ./lib/ncl-2.1.18/ncl/nxsallocatematrix.h ./lib/ncl-2.1.18/ncl/nxscxxdiscretematrix.h ./lib/ncl-2.1.18/ncl/nxsdistancesblock.h ./lib/ncl-2.1.18/ncl/nxssetreader.h ./lib/ncl-2.1.18/ncl/nxstreesblock.h ./lib/ncl-2.1.18/ncl/nxsassumptionsblock.h ./lib/ncl-2.1.18/ncl/nxsdatablock.h ./lib/ncl-2.1.18/ncl/nxsexception.h ./lib/ncl-2.1.18/ncl/nxsstring.h ./lib/ncl-2.1.18/ncl/nxsunalignedblock.h ./lib/ncl-2.1.18/ncl/nxsblock.h ./lib/ncl-2.1.18/ncl/nxsdefs.h ./lib/ncl-2.1.18/ncl/nxsmultiformat.h ./lib/ncl-2.1.18/ncl/nxstaxaassociationblock.h ./lib/ncl-2.1.18/ncl/nxsutilcopy.h ./lib/ncl-2.1.18/ncl/nxscdiscretematrix.h ./lib/ncl-2.1.18/ncl/nxsdiscretedatum.h ./lib/ncl-2.1.18/ncl/nxspublicblocks.h ./lib/ncl-2.1.18/ncl/nxstaxablock.h
libncl_a_CXXFLAGS= -w -O3 -I$(top_srcdir)/lib/ncl-2.1.18/  -D__STRICT_ANSI__ $(AM_CPP_FLAGS)

################################
# cleanup 
# clean : 
# 	-rm -rf $(bin_PROGRAMS) *.o *.a 

doc : 
	doxygen 


man:  
	$(top_srcdir)/manual/process.sh $(top_srcdir)/manual


date :
	$(top_srcdir)/utils/updateDate.sh

mydist : 
	make date 
	make dist 
	make dist-zip
