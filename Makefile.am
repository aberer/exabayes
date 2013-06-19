AUTOMAKE_OPTIONS=foreign

EXTRA_DIST=lib/Random123-1.07 src/pll/avxLikelihood-pll.c src/examl/avxLikelihood.c src/pll/axml-pll.h src/examl/axml-examl.h doc/*

.PHONY: my-update-src doc 

WARN_FLAGS= -Wall --pedantic -Wno-vla -Wno-unused-function -Wno-write-strings  -Wsign-compare
noinst_LIBRARIES=libbayes.a  libncl.a
RAX_FEATURE= -D_GNU_SOURCE -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3 -D_USE_ALLREDUCE -fomit-frame-pointer -funroll-loops

# TODO 
# SIMD_FLAGS is potentially problematic: auto-avx vectorization was considered harmfull? 

# TODO threads
  #  $(srcppre)/genericParallelization-pll.c

bin_PROGRAMS= exabayes prsf asdsf

BAYES_INCLUDES=-Ilib/Random123-1.07/include -Ilib/ncl-2.1.18/ -Isrc/proposals

include src-helper.am
include src.am 

PLLINCLUDE=-I./src/pll
libpll_a_SOURCES= ./src/pll/evaluateGenericSpecial-pll.c ./src/pll/axml-pll.c ./src/pll/trash-pll.c ./src/pll/vec_unit.h ./src/pll/makenewzGenericSpecial-pll.c ./src/pll/mem_alloc-pll.c ./src/pll/restartHashTable-pll.c ./src/pll/axml-pll.h ./src/pll/mem_alloc.h ./src/pll/recom-pll.c ./src/pll/searchAlgo-pll.c ./src/pll/phylip_parser/msa_sites.h ./src/pll/phylip_parser/xalloc.h ./src/pll/phylip_parser/phylip.h ./src/pll/phylip_parser/lexer.h ./src/pll/phylip_parser/ssort.h ./src/pll/utils-pll.c ./src/pll/treeIO-pll.c ./src/pll/cycle.h ./src/pll/aligned_buffer.h ./src/pll/globalVariables.h ./src/pll/topologies-pll.c ./src/pll/utils.h ./src/pll/evaluatePartialGenericSpecial-pll.c ./src/pll/fastDNAparsimony-pll.c ./src/pll/randomTree-pll.c ./src/pll/optimizeModel-pll.c ./src/pll/treeIO.h ./src/pll/newviewGenericSpecial-pll.c ./src/pll/bipartitionList-pll.c ./src/pll/models-pll.c ./src/pll/genericParallelization.h ./src/avxConditional.c
libpll_a_CFLAGS= -w $(SIMD_FLAGS) $(PLLINCLUDE) -O2 $(RAX_FEATURE) 
EXAMLINCLUDE=-I./src/examl/
libexaml_a_SOURCES= ./src/examl/searchAlgo.c ./src/examl/treeIO.c ./src/examl/restartHashTable.c ./src/examl/optimizeModel.c ./src/examl/newviewGenericSpecial.c ./src/examl/trash.c ./src/examl/topologies.c ./src/examl/axml-examl.h ./src/examl/globalVariables.h ./src/examl/evaluatePartialGenericSpecial.c ./src/examl/bipartitionList.c ./src/examl/models.c ./src/examl/makenewzGenericSpecial.c ./src/examl/evaluateGenericSpecial.c ./src/examl/axml.c ./src/avxConditional.c ./src/examl/fastDNAparsimony.c
libexaml_a_CFLAGS= -w $(SIMD_FLAGS)  -fgnu89-inline  $(EXAMLINCLUDE) -O2 $(RAX_FEATURE) 

# configure examl/pll specific stuff 
if WITH_EXAML 
LIBTOUSE=libexaml.a
LIBINCLUDE=$(EXAMLINCLUDE)
noinst_LIBRARIES+=libexaml.a
bin_PROGRAMS+= parserExaml
else
LIBTOUSE=libpll.a
LIBINCLUDE=$(PLLINCLUDE) 
bin_PROGRAMS+= parserPll
noinst_LIBRARIES+=libpll.a
endif 

# AM_CXXFLAGS=-O3
libbayes_a_CXXFLAGS=-O3 $(BAYES_INCLUDES) $(LIBINCLUDE)
exabayes_LDADD= libncl.a libbayes.a 
exabayes_SOURCES = ./src/exabayes.cpp ./src/CommandLine.cpp ./src/CommandLine.hpp 
exabayes_LDADD+=$(LIBTOUSE)
exabayes_CXXFLAGS= -O3 $(BAYES_INCLUDES) $(LIBINCLUDE)

################################################################
# PRSF

prsf_LDADD = libbayes.a $(LIBTOUSE) 
prsf_CXXFLAGS=-O3 $(LIBINCLUDE) # $(BAYES_INCLUDES) 
prsf_SOURCES= ./src/prsf.cpp


################################################################
# ASDSF 
asdsf_CXXFLAGS=-O3  $(LIBINCLUDE) $(BAYES_INCLUDES)
asdsf_LDADD=    libbayes.a $(LIBTOUSE)  libncl.a
asdsf_SOURCES= ./src/asdsf.cpp

################################################################
# examl parser 
pre1=./lib/ExaML/parser/
parserExaml_SOURCES=$(pre1)/axml.c $(pre1)/parsePartitions.c $(pre1)/axml.h  $(pre1)/globalVariables.h



################################################################
# ncl 
npre=lib/ncl-2.1.18/ncl
libncl_a_SOURCES = $(npre)/nxsreader.cpp $(npre)/nxscxxdiscretematrix.cpp $(npre)/nxstaxaassociationblock.cpp $(npre)/nxsunalignedblock.cpp $(npre)/nxsmultiformat.cpp $(npre)/nxsdatablock.cpp $(npre)/nxspublicblocks.cpp $(npre)/nxstaxablock.cpp $(npre)/nxscharactersblock.cpp $(npre)/nxsstring.cpp $(npre)/nxstoken.cpp $(npre)/nxsassumptionsblock.cpp $(npre)/nxsblock.cpp $(npre)/nxsdistancesblock.cpp $(npre)/nxssetreader.cpp $(npre)/nxsexception.cpp $(npre)/nxstreesblock.cpp $(npre)/ncl.h $(npre)/nxscharactersblock.h $(npre)/nxsdistancedatum.h $(npre)/nxsreader.h $(npre)/nxstoken.h $(npre)/nxsallocatematrix.h $(npre)/nxscxxdiscretematrix.h $(npre)/nxsdistancesblock.h $(npre)/nxssetreader.h $(npre)/nxstreesblock.h $(npre)/nxsassumptionsblock.h $(npre)/nxsdatablock.h $(npre)/nxsexception.h $(npre)/nxsstring.h $(npre)/nxsunalignedblock.h $(npre)/nxsblock.h $(npre)/nxsdefs.h $(npre)/nxsmultiformat.h $(npre)/nxstaxaassociationblock.h $(npre)/nxsutilcopy.h $(npre)/nxscdiscretematrix.h $(npre)/nxsdiscretedatum.h $(npre)/nxspublicblocks.h $(npre)/nxstaxablock.h

# noinst_LIBRARIES+=
libncl_a_CXXFLAGS= -w -O3 #  $(AM_CFLAGS) 


################################################################
# pll parser 
pre2=./lib/phylogenetic-likelihood-library/parser
parserPll_SOURCES=$(pre2)/axml.c $(pre2)/parsePartitions.c $(pre2)/axml.h $(pre2)/globalVariables.h src/pll/mem_alloc.h ./src/pll/mem_alloc-pll.c
parserPll_CFLAGS=  -w -I$(pre2) $(PLLINCLUDE) $(AM_CFLAGS) 

################################
# cleanup 

clean:
	-rm -rf  $(bin_PROGRAMS) $(exabayes_OBJECTS) $(libbayes_a_OBJECTS) libbayes.a

clean-all : 
	-rm -rf $(bin_PROGRAMS) *.o *.a 


doc : 
	doxygen 
