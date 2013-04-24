/**
   @file bayes.c
 
   @brief The top-level functionality of ExaBayes.

*/


#include <sstream>

#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"
#include "proposals.h"
#include "nclConfigReader.h"
#include "misc-utils.h"
#include "chain.h"
#include "adapters.h"
#include "eval.h"
#include "proposals.h"
#include "tune.h"
#include "prsfComputer.h"
#include "TreeAln.hpp"
#include "AvgSplitFreqAssessor.hpp"

extern double masterTime; 



/** 
    @brief nasty, but we need to manage a lot of global
    variables. 

    Maybe we can turn this into a singleton later or
    de-globalize it. Has already caused some trouble. 

 */ 
void setupGlobals(initParamStruct *initParams)
{
  if (initParams->numGen > 0)
    gAInfo.numGen = initParams->numGen; 

  gAInfo.samplingFrequency = initParams->samplingFrequency; 
  gAInfo.diagFreq = initParams->diagFreq; 
  gAInfo.numberOfRuns =   initParams->numIndiChains; 
  gAInfo.numberCoupledChains = initParams->numCoupledChains; 

  gAInfo.printFreq = initParams->printFreq; 
  gAInfo.asdsfIgnoreFreq = initParams->asdsfIgnoreFreq; 
  gAInfo.asdsfConvergence = initParams->asdsfConvergence; 
  gAInfo.heatFactor  = initParams->heatFactor; 
  gAInfo.swapInterval =  initParams->swapInterval; 
  gAInfo.tuneHeat = initParams->tuneHeat; 
  gAInfo.burninGen = initParams->burninGen; 
  gAInfo.burninProportion = initParams->burninProportion; 
  gAInfo.tuneFreq = initParams->tuneFreq ; 

  /* initialize a matrix of swaps (wasting some space here) */
  // gAInfo.swapInfo = (SuccessCtr**)exa_calloc(gAInfo.numberOfRuns, sizeof(SuccessCtr*)); 
  // int n = gAInfo.numberCoupledChains; 
  // for(int i = 0; i < gAInfo.numberOfRuns; ++i)
  //   gAInfo.swapInfo[i] = (SuccessCtr*)exa_calloc( n * n , sizeof(SuccessCtr)); 

  // gAInfo.temperature = (double*)exa_calloc(gAInfo.numberOfRuns, sizeof(double)); 
  // for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    // gAInfo.temperature[i] = gAInfo.heatFactor; 
}




#include <vector>

#include "CoupledChains.hpp"

//STAY 
bool convergenceDiagnostic(vector<CoupledChains*> runs)
{
  Chain *allChains = NULL; 
  assert(0); 
  // just force to compile 


  if(gAInfo.numberOfRuns > 1)
    { 
      vector<string> fns; 
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	{
	  stringstream ss; 
	  ss <<  PROGRAM_NAME << "_topologies." << run_id << "." << i; 
	  fns.push_back(ss.str());
	}
     
      AvgSplitFreqAssessor asdsf(fns);

      int end = asdsf.getEnd();
      
      int treesInBatch = gAInfo.diagFreq / gAInfo.samplingFrequency; 

      end /= treesInBatch; 
      end *= treesInBatch;       

      if(end > 0)
	{	  
	  asdsf.setEnd(end);
	  if( gAInfo.burninGen > 0 )
	    {
	      assert(gAInfo.burninProportion == 0.); 

	      int treesToDiscard =  gAInfo.burninGen / gAInfo.samplingFrequency; 

	      if(end < treesToDiscard + 2 )
		return false; 
	      else 
		asdsf.setStart(treesToDiscard);  
	    }
	  else 
	    {
	      assert(gAInfo.burninGen == 0); 
	      int start = (int)((double)end * gAInfo.burninProportion  ); 
	      asdsf.setStart(start);
	    } 

	  asdsf.extractBips();
	  double asdsfVal = asdsf.computeAsdsf(gAInfo.asdsfIgnoreFreq);

	  PRINT("ASDSF for trees %d-%d: %f\n", asdsf.getStart(), asdsf.getEnd(), asdsfVal ); 

	  return asdsfVal < gAInfo.asdsfConvergence; 

	}
      else 
	return false; 
      
    }
  else 
    return allChains[0].currentGeneration > gAInfo.numGen; 
}



/**
   @brief run all chains for a given number of generations
 */
void runChains(vector<CoupledChains*> allRuns, int diagFreq)
{  
  bool hasConverged = false;   
  while(not hasConverged)
    {      
      for(nat i = 0; i < allRuns.size(); ++i)
	{
	  auto run = allRuns[i]; 
	  run->executePart(diagFreq);
	}

      hasConverged = convergenceDiagnostic(allRuns); 

#ifdef ENABLE_PRSF
      if(isOutputProcess)
	printPRSF(run_id);
#endif
    }
}



// #define TEST 
 
#include "MyTestProposal.hpp"



/**
   @brief the main ExaBayes function.

   @param tr -- a tree structure that has been initialize in one of the adapter mains. 
   @param adef -- the legacy adef
 */
void exa_main (analdef *adef, int seed, initParamStruct *initParams)
{   
  Chain *indiChains = NULL; 		/* one state per indipendent run/chain */  

  timeIncrement = gettime();
  gAInfo.adef = adef; 

#ifdef TEST   
  MyTestProposal prop();
  
  

  exit(0); 
#endif

  vector<CoupledChains*> runs; 

  initializeIndependentChains(adef,  seed, runs, initParams); 
  assert(gAInfo.numberCoupledChains > 0);

  assert(gAInfo.diagFreq != 0 ); 
  runChains(runs, gAInfo.diagFreq); 

  if(isOutputProcess() )
    {
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	finalizeOutputFiles(indiChains + i);
      PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    } 
}

