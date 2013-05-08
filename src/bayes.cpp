/**
   @file bayes.c
 
   @brief The top-level functionality of ExaBayes.

*/


#include <sstream>
#include <vector>

#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "output.h"
#include "proposals.h"
#include "nclConfigReader.h"
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

}


//STAY 
bool convergenceDiagnostic(vector<CoupledChains*> runs)
{
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

 #if HAVE_PLL == 0      
	  if(processID == 0)
#endif 
	    cout << "ASDSF for trees " << asdsf.getStart() << "-" << asdsf.getEnd() << ": " <<setprecision(2) << asdsfVal * 100 << "%" << endl; 

	  return asdsfVal < gAInfo.asdsfConvergence; 

	}
      else 
	return false; 
      
    }
  else 
    {
      return runs[0]->getChain(0)->currentGeneration > gAInfo.numGen;
    }
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


// LEGACY
static void traverseInitFixedBL(nodeptr p, int *count, TreeAln *traln,  double z )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;
  
  for( i = 0; i < traln->getNumBranches(); i++)
      traln->setBranchLengthSave(z, i, p); 
  
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitFixedBL(q->back, count, traln, z);
	  q = q->next;
	} 
    }
}


static void initTreeWithOneRandom(int seed, vector<TreeAln*> &tralns)
{

  TreeAln *traln = tralns[0]; 
  tree *tr = traln->getTr();

  TreeRandomizer trRandomizer(seed, traln);
  trRandomizer.randomizeTree();

  int count = 0; 
  traverseInitFixedBL( tr->start->back, &count, traln, TreeAln::initBL);
  assert(count  == 2 * tr->mxtips - 3);
  
  for(nat i = 1; i < tralns.size(); ++i)
    *(tralns[i])  = *traln;   
}


static void initWithStartingTree(FILE *fh, vector<TreeAln*> &tralns)
{
  // fetch a tree 
  TreeAln *traln = tralns[0]; 
  tree *tr = traln->getTr();
  boolean hasBranchLength =  readTreeWithOrWithoutBL(tr, fh);

  int count = 0;       
  if(hasBranchLength)
    traverseInitCorrect(tr->start->back, &count, traln ) ;  
  else      
    traverseInitFixedBL(tr->start->back, &count, traln, TreeAln::initBL ); 

  assert(count == 2 * tr->mxtips  -3);       

  for(nat i = 1 ; i < tralns.size(); ++i)
    *(tralns[i]) = *traln; 
}


static int countNumberOfTreesQuick(char *fn )
{
  FILE *fh = fopen(fn, "r"); 

  if(fh == 0)
    return 0; 

  int c = 0; 
  int result = 0; 
  while( ( c = getc(fh) ) != EOF)
    {
      if(c == ';')
	++result; 	
    }
  
  fclose(fh);
  return result; 
}


static void initializeIndependentChains( analdef *adef, int seed, vector<CoupledChains*> &runs, initParamStruct *initParams )
{
  FILE *treeFH = NULL; 
  int numTrees = countNumberOfTreesQuick(tree_file); 
  
  if( numTrees > 0 )
    treeFH = myfopen(tree_file, "r"); 

  PRINT("number of independent runs=%d, number of coupled chains per run=%d \n", gAInfo.numberOfRuns, gAInfo.numberCoupledChains ); 

#ifdef DEBUG_LNL_VERIFY
  gAInfo.debugTree = new TreeAln();   
  gAInfo.debugTree->initializeFromByteFile(byteFileName); 
  gAInfo.debugTree->enableParsimony();
#endif

  // sets up tree structures 
  vector<TreeAln*>  trees; 
  for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
    {
      TreeAln *traln = new TreeAln();
      traln->initializeFromByteFile(byteFileName); 
#if HAVE_PLL != 0 
      traln->enableParsimony();
#endif
      trees.push_back(traln); 
    }


  Randomness masterRand(seed);   
  vector<int> runSeeds; 
  vector<int> treeSeeds; 
  for(int i = 0; i < gAInfo.numberOfRuns;++i)
    {
      randCtr_t r = masterRand.generateSeed();
      runSeeds.push_back(r.v[0]); 
      treeSeeds.push_back(r.v[1]); 
    }


  for(int i = 0; i < gAInfo.numberOfRuns ; ++i)
    {      
      if( i < numTrees)
	initWithStartingTree(treeFH, trees); 
      else 
	initTreeWithOneRandom(treeSeeds[i], trees);

#if HAVE_PLL == 0
      if(i % initParams->numRunParallel != gAInfo.myBatch )
	continue; 
#endif

      CoupledChains *cc = new CoupledChains(runSeeds[i], gAInfo.numberCoupledChains, trees, i, initParams);
      runs.push_back(cc); 
    }
  
  // only use one restorer for all chains 
  LnlRestorer *restorer = new LnlRestorer(runs[0]->getChain(0));
  for(auto run : runs )
    {
      for(int i = 0; i < run->getNumberOfChains(); ++i )
	{
	  Chain *chain = run->getChain(i); 
	  chain->setRestorer(restorer); 
	}
    }

  if(numTrees > 0)
    fclose(treeFH); 

  PRINT("\n"); 

}




// #include "branch.h"
// #include "Path.hpp"
#include "TreeRandomizer.hpp" 



// #define TEST  

// #if
// extern "C"
// {
//   unsigned int evaluateParsimony(tree *tr, partitionList *pr, nodeptr p, boolean full); 
  
// }

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

  vector<CoupledChains*> runs; 

  initializeIndependentChains(adef,  seed, runs, initParams); 
  assert(gAInfo.numberCoupledChains > 0);

  // Chain* chain = runs[0]->getChain(0); 
  // cout << setprecision(10) << "interanl=" <<  TreeAln::zMax << "\t" <<  TreeAln::zMin << endl; 
  // cout << "min=" <<  branchLengthToReal(chain->traln->getTr(), TreeAln::zMin) << 
  //   "\tmax=" << branchLengthToReal(chain->traln->getTr(), TreeAln::zMax) << endl;  
  // assert(0);
  

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

