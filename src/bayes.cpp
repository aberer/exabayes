/**
   @file bayes.c
 
   @brief The top-level functionality of ExaBayes.

*/


#include <sstream>
#include <vector>

#include "axml.h"
#include "bayes.h"
#include "GlobalVariables.hpp"
#include "output.h"
#include "proposals.h"
#include "ConfigReader.hpp"
#include "misc-utils.h"
#include "adapters.h"
#include "eval.h"
#include "proposals.h"
#include "prsfComputer.h"
#include "TreeAln.hpp"
#include "AvgSplitFreqAssessor.hpp"
#include "LnlRestorer.hpp"
#include "TreeRandomizer.hpp"
#include "treeRead.h"
#include "CoupledChains.hpp"
#include "SampleMaster.hpp"

extern double masterTime; 



/** 
    @brief nasty, but we need to manage a lot of global
    variables. 

    Maybe we can turn this into a singleton later or
    de-globalize it. Has already caused some trouble. 

 */ 
void setupGlobals(initParamStruct *initParams)
{
  // if (initParams->numGen > 0)
    // globals.numGen = initParams->numGen; 

  // globals.samplingFrequency = initParams->samplingFrequency; 
  // globals.diagFreq = initParams->diagFreq; 
  // globals.numberOfRuns =   initParams->numIndiChains; 
  // globals.numberCoupledChains = initParams->numCoupledChains; 
  // globals.asdsfIgnoreFreq = initParams->asdsfIgnoreFreq; 
  // globals.asdsfConvergence = initParams->asdsfConvergence; 
  // globals.burninGen = initParams->burninGen; 
  // globals.burninProportion = initParams->burninProportion; 
}

// #define TEST  



/**
   @brief the main ExaBayes function.

  @param tr -- a tree structure that has been initialize in one of the adapter mains. 
   @param adef -- the legacy adef
 */
void exa_main (analdef *adef, int seed, initParamStruct *initParams)
{   
  timeIncrement = gettime();
  globals.adef = adef; 

#ifdef TEST   
  TreeAln traln; 
  traln.initializeFromByteFile(byteFileName);
  traln.enableParsimony(); 

  TreeRandomizer r(123, &traln); 
  r.randomizeTree();
  tree *tr = traln.getTr(); 
  
  for(int i = tr->mxtips+1; i < 2 * tr->mxtips; ++i)
    {
      cout << exa_evaluateParsimony(traln, tr->nodep[i], TRUE ) << endl; 
    }

  exit(0);

#endif


  SampleMaster master(adef,  seed, initParams);

  master.run();


  // if(isOutputProcess() )
  //   {
  //     for(nat i = 0; i < runs.size(); ++i)
  // 	finalizeOutputFiles(indiChains + i);
  //     PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
  //     PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
  //   }
  
}

