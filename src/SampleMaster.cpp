#include <sstream>

#include "output.h"
#include "eval.h"
#include "adapters.h"
#include "SampleMaster.hpp"
#include "Chain.hpp"
#include "TreeRandomizer.hpp"
#include "treeRead.h"

// TODO =( 
#include "GlobalVariables.hpp"

#include "LnlRestorer.hpp"
#include "AvgSplitFreqAssessor.hpp"



extern double masterTime; 

//STAY 
bool SampleMaster::convergenceDiagnostic()
{
  if(numRunConv > 1)    
    { 
      vector<string> fns; 
      for(int i = 0; i < numRunConv; ++i)
	{
	  stringstream ss; 
	  ss <<  PROGRAM_NAME << "_topologies." << runId << "." << i; 
	  fns.push_back(ss.str());
	}
     
      AvgSplitFreqAssessor asdsf(fns);

      int end = asdsf.getEnd();
      
      int treesInBatch = diagFreq / samplingFreq; 

      end /= treesInBatch; 
      end *= treesInBatch;       

      if(end > 0)
	{	  
	  asdsf.setEnd(end);
	  if( burninGen > 0 )
	    {
	      assert(burninProportion == 0.); 

	      int treesToDiscard =  burninGen / samplingFreq; 

	      if(end < treesToDiscard + 2 )
		return false; 
	      else 
		asdsf.setStart(treesToDiscard);  
	    }
	  else 
	    {
	      assert(burninGen == 0); 
	      int start = (int)((double)end * burninProportion  ); 
	      asdsf.setStart(start);
	    } 

	  asdsf.extractBips();
	  double asdsfVal = asdsf.computeAsdsf(asdsfIgnoreFreq);

#if HAVE_PLL == 0      
	  if(processID == 0)
#endif 
	    cout << "ASDSF or trees " << asdsf.getStart() << "-" << asdsf.getEnd() << ": " <<setprecision(2) << asdsfVal * 100 << "%" << endl; 

	  return asdsfVal < asdsfConvergence; 

	}
      else 
	return false; 
      
    }
  else 
    {
      return runs[0].getChain(0)->currentGeneration > numGen;
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


// TODO finalize output files  

static int countNumberOfTreesQuick(const char *fn )
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


// void SampleMaster::initializeRunParameters()
// {
  

// }





SampleMaster::SampleMaster(const CommandLine &cl , ParallelSetup &pl) 
  : 
  runId(cl.getRunid())
  // : diagFreq(initParams->diagFreq)
  // , asdsfIgnoreFreq(initParams->asdsfIgnoreFreq)
  // , asdsfConvergence(initParams->asdsfConvergence)
  // , burninGen(initParams->burninGen)
  // , burninProportion(initParams->burninProportion)
  // , samplingFreq(initParams->samplingFrequency)
  // , numRunConv(initParams->numIndiChains)
  // , numGen(initParams->numGen)
  // , myBatch(_myBatch)
{
  ConfigReader reader; 
  
  initParamStruct *initParams; 
  assert(0); 
  
  assert(initParams->numCoupledChains != 0 ); 

  FILE *treeFH = NULL; 
  int numTrees = countNumberOfTreesQuick(cl.getTreeFile().c_str()); 
  
  if( numTrees > 0 )
    treeFH = myfopen(cl.getTreeFile().c_str(), "r"); 

  PRINT("number of independent runs=%d, number of coupled chains per run=%d \n", initParams->numIndiChains, initParams->numCoupledChains ); 

#ifdef DEBUG_LNL_VERIFY
  globals.debugTree = new TreeAln();   
  globals.debugTree->initializeFromByteFile(byteFileName); 
  globals.debugTree->enableParsimony();
#endif

  // sets up tree structures 
  vector<TreeAln*>  trees; 
  
  for(int i = 0; i < initParams->numCoupledChains; ++i)
    {
      TreeAln *traln = new TreeAln();
      traln->initializeFromByteFile(cl.getAlnFileName()); 
#if HAVE_PLL != 0 
      traln->enableParsimony();
#endif
      trees.push_back(traln); 
    }


  Randomness masterRand(cl.getSeed());   
  vector<int> runSeeds; 
  vector<int> treeSeeds; 
  for(int i = 0; i < initParams->numIndiChains;++i)
    {
      randCtr_t r = masterRand.generateSeed();
      runSeeds.push_back(r.v[0]); 
      treeSeeds.push_back(r.v[1]); 
    }


  for(int i = 0; i < initParams->numIndiChains ; ++i)
    {      
      if( i < numTrees)
	initWithStartingTree(treeFH, trees); 
      else 
	initTreeWithOneRandom(treeSeeds[i], trees);

#if HAVE_PLL == 0
      if(i % initParams->numRunParallel != myBatch )
	continue; 
#endif

      runs.push_back(CoupledChains(runSeeds[i], initParams->numCoupledChains, trees, i, initParams));       
    }
  
  if(initParams->tuneHeat)
    for(auto r : runs)
      r.enableHeatTuning(initParams->tuneFreq); 
  
  // only use one restorer for all chains 
  LnlRestorer *restorer = new LnlRestorer(runs[0].getChain(0));
  for(auto run : runs )
    {
      for(int i = 0; i < run.getNumberOfChains(); ++i )
	{
	  Chain *chain = run.getChain(i); 
	  chain->setRestorer(restorer); 
	}
    }

  if(numTrees > 0)
    fclose(treeFH); 

  PRINT("\n"); 

}


void SampleMaster::run()
{
  bool hasConverged = false;   
  while(not hasConverged)   
    {      
      for(nat i = 0; i < runs.size(); ++i)
	{
	  auto run = runs[i]; 
	  run.executePart(diagFreq);
	}

      hasConverged = convergenceDiagnostic(); 

#ifdef ENABLE_PRSF
      if(isOutputProcess())
	printPRSF(run_id);
#endif
    }
}
 
 
void SampleMaster::finalizeRuns()
{
  if(isOutputProcess() )
    {
      for(auto run : runs)
	{
	  for(int i = 0; i < run.getNumberOfChains(); ++i)
	    {
	      Chain *chain = run.getChain(0);
	      if( chain->getChainHeat() == 1.f)
		finalizeOutputFiles(chain);
	    }
	}

      PRINT("\nConverged after %d generations\n",  runs[0].getChain(0)->getGeneration());
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    }  
}
