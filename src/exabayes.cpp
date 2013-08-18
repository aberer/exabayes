/** 
    @file exabayes.cpp
    
    @brief This file sets the flavour of exabayes that has been
    compiled (i.e., sequential, pll, examl, ...)
    
*/ 



// TODO re-activate all that initial bla bla when starting up
// exa-bayes (model / program info )


#ifdef HAVE_AVX
#define __AVX
#endif

#include <iostream>
#include <sstream>
#include <chrono>

#include "axml.h" 

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "time.hpp"

#include "config/CommandLine.hpp"
#include "SampleMaster.hpp"
#include "ParallelSetup.hpp"

#include "teestream.hpp"

#define TEST  

#ifdef TEST
#include "parameters/BranchLengthsParameter.hpp"
#include "TreeRandomizer.hpp"
#include "Chain.hpp"
#include "BranchLengthMultiplier.hpp"
#include "ParsimonyEvaluator.hpp"
#include "BipartitionHashNew.hpp"
#endif


// have ae look at that later again 
double fastPow(double a, double b)
{
  union 
  {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}


/**
   @brief the main ExaBayes function.

  @param tr -- a tree structure that has been initialize in one of the adapter mains. 
   @param adef -- the legacy adef
 */
static void exa_main (const CommandLine &cl, const ParallelSetup &pl )
{   
  timeIncrement = CLOCK::system_clock::now(); 

#ifdef TEST     
  // auto traln = TreeAln{};
  TreeAln traln; 
  traln.initializeFromByteFile(cl.getAlnFileName());
  auto seed = randCtr_t{123,456};
  Randomness rand(seed); 

  auto map = traln.getNameMap(); 
  std::cout << "map: " << map << std::endl; 

  TreeRandomizer::randomizeTree(traln, rand); 
  auto blParam = BranchLengthsParameter(0,0);
  blParam.addPartition(0);
  auto ap = dynamic_cast<AbstractParameter*>(&blParam); 
  auto branches = traln.extractBranches(ap); 

  double ctr = 0.001; 
  for(auto &b : branches)
    {
      b.setLength(ctr, ap);
      traln.setBranch(b, ap);
      ctr += 0.001;
    }

  auto theHash = BipartitionHashNew{traln.getNumberOfTaxa()}; 

  theHash.addTree(traln, true);
  TreeRandomizer::randomizeTree(traln,rand);
  theHash.addTree(traln,true);
  theHash.addTree(traln,true);

  for(auto &elem : theHash)
    {
      auto &bip = elem.first; 
      std::cout << bip << "\t" << elem.second << std::endl; 
      std::cout << theHash.getBranchLengths(elem.first) << std::endl;
      std::cout << "\t=>" << bip.count() << std::endl; 
      std::cout << "full: "; 
      bip.printVerbose(std::cout, map);
      std::cout << std::endl; 
    }

  exit(0); 
#else 
  SampleMaster master(  pl, cl );
  master.initializeRuns(); 
  master.run();
  master.finalizeRuns();
#endif
}


/* 
   tell the CPU to ignore exceptions generated by denormalized floating point values.
   If this is not done, depending on the input data, the likelihood functions can exhibit 
   substantial run-time differences for vectors of equal length.
*/

void ignoreExceptionsDenormFloat()
{
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
  _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif   
}


#if HAVE_PLL != 0
#include "globalVariables.h"
#endif


void makeInfoFile(const CommandLine &cl, const ParallelSetup &pl )
{
  stringstream ss; 

  ss << OutputFile::getFileBaseName(cl.getWorkdir()) << "_info."  << cl.getRunid() ;

  if( pl.isGlobalMaster() && std::ifstream(ss.str()))
    {
      std::cerr << pl << std::endl; 
      std::cerr << std::endl <<  "File " << ss.str() << " already exists (probably \n"
		<< "from previous run). Please choose a new run-id or remove previous output files. " << std::endl; 
      ParallelSetup::genericExit(-1); 
    }

  globals.logFile = ss.str();   
  
  globals.logStream =  new ofstream  (globals.logFile); 
  globals.teeOut =  new teestream(cout, *globals.logStream);

  if(not pl.isGlobalMaster())
    tout.disable(); 
}


void initializeProfiler()
{
  // see this page for info 
  // http://google-perftools.googlecode.com/svn/trunk/doc/cpuprofile.html  
  // that option is important
  // CPUPROFILE_FREQUENCY=x
#ifdef _USE_GOOGLE_PROFILER
  ProfilerStart("profile.out");
#endif
}
 

 
void finalizeProfiler()
{
#ifdef _USE_GOOGLE_PROFILER
  ProfilerStop();
#endif
}


int main(int argc, char **argv)
{ 
  ParallelSetup pl(argc,argv); 		// MUST be the first thing to do because of mpi_init ! 

  initializeProfiler();

#if HAVE_PLL != 0 && ( (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS)))
  assert(0); 
#endif

  ignoreExceptionsDenormFloat(); 
  CommandLine cl(argc, argv); 

#if HAVE_PLL == 0 
  pl.initializeExaml(cl);
#endif

  makeInfoFile(cl, pl);

  // cl.printVersion(true);  
  tout << endl; 

  tout << "This is " << PROGRAM_NAME << " version "
       << VERSION
       << ", a tool for Bayesian MCMC sampling of phylogenetic trees." 
#if HAVE_PLL == 0 
       << "\nbuild with the ExaML code base and MPI-support."
#else 
       << "\nbuild with the (Phylogenetic Likelihood Library) PLL code base for sequential execution."
#endif
       << "\nPlease send any feature requests and inquiries to exabayes-at-googlemail-dot-com"    
       << std::endl << std::endl; 

  tout << PROGRAM_NAME << " was called as follows: " << endl; 
  for(int i = 0; i < argc; ++i)
    tout << argv[i] << " " ; 
  tout << endl << endl; 

  exa_main( cl, pl); 

  finalizeProfiler();
  pl.finalize();  
  
  delete globals.logStream; 
  delete globals.teeOut;
  return 0;
}
