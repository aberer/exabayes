AUTOMAKE_OPTIONS=foreign

EXTRA_DIST=lib/Random123-1.08 src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c src/pll/axml-pll.h src/examl/axml-examl.h examples/* build.sh apple-build.sh $(top_srcdir)/manual/manual.pdf manual/css/org.css manual/img/banner.png $(top_srcdir)/manual/manual.html README.txt LICENSE.txt COPYING

.PHONY: my-update-src doc mydist  man  dist

WARN_FLAGS= -Wall --pedantic -Wextra -Wno-unused-function -Wno-unused-parameter  -Wuninitialized -Wno-mismatched-tags -Wno-unused-private-field
noinst_LIBRARIES=libbayes.a libncl.a 
RAX_FEATURE= -D_GNU_SOURCE -D__SIM_SSE3 -D_USE_ALLREDUCE -fomit-frame-pointer -funroll-loops -D__STRICT_ANSI__ -D_OPTIMIZED_FUNCTIONS

MY_SIMD_FLAGS= $(SIMD_FLAGS)
# MY_SIMD_FLAGS= -msse4.2 -DMANUAL_AVX_OVERRIDE

AM_CPP_FLAGS=-fno-common

SRCINCLUDE=-I$(top_srcdir)/src

BAYES_INCLUDES=-I$(top_srcdir)/lib/Random123-1.08/include -I$(top_srcdir)/lib/ncl-2.1.18/  $(SRCINCLUDE) 

include src-helper.am
include src.am 

PLLDIR=$(top_srcdir)/src/pll
libpll_a_SOURCES= $(PLLDIR)/evaluateGenericSpecial-pll.c $(PLLDIR)/axml-pll.c $(PLLDIR)/trash-pll.c $(PLLDIR)/vec_unit.h $(PLLDIR)/makenewzGenericSpecial-pll.c $(PLLDIR)/mem_alloc-pll.c $(PLLDIR)/restartHashTable-pll.c $(PLLDIR)/axml-pll.h $(PLLDIR)/mem_alloc.h $(PLLDIR)/recom-pll.c $(PLLDIR)/searchAlgo-pll.c $(PLLDIR)/phylip_parser/msa_sites.h $(PLLDIR)/phylip_parser/xalloc.h $(PLLDIR)/phylip_parser/phylip.h $(PLLDIR)/phylip_parser/lexer.h $(PLLDIR)/phylip_parser/ssort.h $(PLLDIR)/utils-pll.c $(PLLDIR)/treeIO-pll.c $(PLLDIR)/cycle.h $(PLLDIR)/aligned_buffer.h $(PLLDIR)/globalVariables.h $(PLLDIR)/topologies-pll.c $(PLLDIR)/utils.h $(PLLDIR)/evaluatePartialGenericSpecial-pll.c $(PLLDIR)/fastDNAparsimony-pll.c $(PLLDIR)/randomTree-pll.c $(PLLDIR)/optimizeModel-pll.c $(PLLDIR)/treeIO.h $(PLLDIR)/newviewGenericSpecial-pll.c $(PLLDIR)/bipartitionList-pll.c $(PLLDIR)/models-pll.c $(PLLDIR)/genericParallelization.h $(top_srcdir)/src/avxConditional.c
libpll_a_CFLAGS= -w $(MY_SIMD_FLAGS) $(RAX_FEATURE) -O2 $(SRCINCLUDE)  $(AM_CPP_FLAGS)
EXAMLDIR=$(top_srcdir)/src/examl
libexaml_a_SOURCES= $(EXAMLDIR)/searchAlgo.c $(EXAMLDIR)/treeIO.c $(EXAMLDIR)/restartHashTable.c $(EXAMLDIR)/optimizeModel.c $(EXAMLDIR)/newviewGenericSpecial.c $(EXAMLDIR)/trash.c $(EXAMLDIR)/topologies.c $(EXAMLDIR)/axml-examl.h $(EXAMLDIR)/globalVariables.h $(EXAMLDIR)/evaluatePartialGenericSpecial.c $(EXAMLDIR)/bipartitionList.c $(EXAMLDIR)/models.c $(EXAMLDIR)/makenewzGenericSpecial.c $(EXAMLDIR)/evaluateGenericSpecial.c $(EXAMLDIR)/axml.c $(top_srcdir)/src/avxConditional.c $(EXAMLDIR)/fastDNAparsimony.c $(EXAMLDIR)/mpiMock.c $(EXAMLDIR)/mpiMock.h
libexaml_a_CFLAGS= -w $(MY_SIMD_FLAGS)  -fgnu89-inline   $(RAX_FEATURE) -O2  $(SRCINCLUDE)  $(AM_CPP_FLAGS)

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

libbayes_a_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(MY_SIMD_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)
libbayes_a_CFLAGS=-O2  $(BAYES_INCLUDES) $(MY_SIMD_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)
exabayes_LDADD= libncl.a libbayes.a 
exabayes_LDADD+=$(LIBTOUSE)
exabayes_SOURCES = $(top_srcdir)/src/apps/exabayes.cpp $(top_srcdir)/src/config/CommandLine.cpp $(top_srcdir)/src/config/CommandLine.hpp
exabayes_CXXFLAGS=-O3  $(BAYES_INCLUDES) $(MY_SIMD_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)


yggdrasil_LDADD= libncl.a libbayes.a 
yggdrasil_LDADD+=$(LIBTOUSE)
yggdrasil_SOURCES = $(top_srcdir)/src/apps/exabayes.cpp $(top_srcdir)/src/config/CommandLine.cpp $(top_srcdir)/src/config/CommandLine.hpp
yggdrasil_CXXFLAGS=-O3  $(BAYES_INCLUDES) $(MY_SIMD_FLAGS) $(RAX_FEATURE) $(WARN_FLAGS) $(AM_CPP_FLAGS)


if WITH_EXAML
cleanobjs=$(exabayes_OBJECTS)
else 
cleanobjs=$(yggdrasil_OBJECTS)
endif


################################################################
# ASDSF 
asdsf_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(MY_SIMD_FLAGS) $(RAX_FEATURE) $(AM_CPP_FLAGS)
asdsf_LDADD= libbayes.a  $(LIBTOUSE)  libncl.a
asdsf_SOURCES= $(top_srcdir)/src/apps/asdsf.cpp 

################################################################
# the parser 
parser_SOURCES=$(top_srcdir)/src/apps/parser.cpp
parser_LDADD= libbayes.a  $(LIBTOUSE) libncl.a
parser_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(MY_SIMD_FLAGS) $(RAX_FEATURE) -Wno-deprecated-writable-strings  $(AM_CPP_FLAGS)

################################################################
# extractBips 
extractBips_CXXFLAGS=-O3 $(MY_SIMD_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
extractBips_LDADD= libbayes.a $(LIBTOUSE)  libncl.a
extractBips_SOURCES=$(top_srcdir)/src/apps/extractBips.cpp  


################################################################
# consense 
consense_CXXFLAGS=-O3 $(MY_SIMD_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
consense_LDADD= libbayes.a  $(LIBTOUSE)  libncl.a
consense_SOURCES=$(top_srcdir)/src/apps/consense.cpp


################################################################
# credibleSet 
credibleSet_CXXFLAGS=-O3 $(MY_SIMD_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
credibleSet_LDADD= libbayes.a  $(LIBTOUSE) libncl.a
credibleSet_SOURCES=$(top_srcdir)/src/apps/credibleSet.cpp

################################################################
# postProcParam 
postProcParam_CXXFLAGS=-O3 $(MY_SIMD_FLAGS) $(BAYES_INCLUDES) $(RAX_FEATURE) $(AM_CPP_FLAGS)
postProcParam_LDADD= libbayes.a  $(LIBTOUSE) 
postProcParam_SOURCES=$(top_srcdir)/src/apps/postProcParam.cpp

################################################################
# ncl 
npre=$(top_srcdir)/lib/ncl-2.1.18/ncl
libncl_a_SOURCES = $(npre)/nxsreader.cpp $(npre)/nxscxxdiscretematrix.cpp $(npre)/nxstaxaassociationblock.cpp $(npre)/nxsunalignedblock.cpp $(npre)/nxsmultiformat.cpp $(npre)/nxsdatablock.cpp $(npre)/nxspublicblocks.cpp $(npre)/nxstaxablock.cpp $(npre)/nxscharactersblock.cpp $(npre)/nxsstring.cpp $(npre)/nxstoken.cpp $(npre)/nxsassumptionsblock.cpp $(npre)/nxsblock.cpp $(npre)/nxsdistancesblock.cpp $(npre)/nxssetreader.cpp $(npre)/nxsexception.cpp $(npre)/nxstreesblock.cpp $(npre)/ncl.h $(npre)/nxscharactersblock.h $(npre)/nxsdistancedatum.h $(npre)/nxsreader.h $(npre)/nxstoken.h $(npre)/nxsallocatematrix.h $(npre)/nxscxxdiscretematrix.h $(npre)/nxsdistancesblock.h $(npre)/nxssetreader.h $(npre)/nxstreesblock.h $(npre)/nxsassumptionsblock.h $(npre)/nxsdatablock.h $(npre)/nxsexception.h $(npre)/nxsstring.h $(npre)/nxsunalignedblock.h $(npre)/nxsblock.h $(npre)/nxsdefs.h $(npre)/nxsmultiformat.h $(npre)/nxstaxaassociationblock.h $(npre)/nxsutilcopy.h $(npre)/nxscdiscretematrix.h $(npre)/nxsdiscretedatum.h $(npre)/nxspublicblocks.h $(npre)/nxstaxablock.h
libncl_a_CXXFLAGS= -w -O3 -I$(top_srcdir)/lib/ncl-2.1.18/  -D__STRICT_ANSI__ $(AM_CPP_FLAGS)

################################
# cleanup 
clean : 
	-rm -rf $(bin_PROGRAMS) *.o *.a 

doc : 
	doxygen 


man:  
	./manual/process.sh


date :
	./utils/updateDate.sh

mydist : 
	make man
	make date 
	make dist 
