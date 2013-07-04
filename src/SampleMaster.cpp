#include <sstream>
#include <fstream>
#include <memory>

#include "config/BlockProposalConfig.hpp"
#include "config/BlockRunParameters.hpp"
#include "config/ConfigReader.hpp"

#include "SampleMaster.hpp"
#include "Chain.hpp"
#include "TreeRandomizer.hpp"
#include "treeRead.h"
#include "tune.h"
#include "RunFactory.hpp"	

#include "GlobalVariables.hpp"

#include "LnlRestorer.hpp"
#include "AvgSplitFreqAssessor.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"

static int countNumberOfTreesQuick(const char *fn ); 
static void initWithStartingTree(FILE *fh, vector<shared_ptr<TreeAln>  > &tralns); 

SampleMaster::SampleMaster(const ParallelSetup &pl, const CommandLine& _cl ) 
  : pl(pl)
  , initTime(CLOCK::system_clock::now())
  , masterRand(cl.getSeed())
  , cl(_cl)
{
  // DO NOT DELETE THIS LINE: 
  // it has cost me ~ 2 days, to figure out that this stupid line
  // throws an increadible number of memory problems
  masterRand = Randomness(cl.getSeed()); 
}


void SampleMaster::initTrees(vector<shared_ptr<TreeAln> > &trees, randCtr_t seed, nat &treesConsumed, nat numTreesAvailable, FILE *treeFH)
{
  Randomness treeRandomness(seed);
  
  vector<shared_ptr<TreeAln> > treesToInitialize; 
  treesToInitialize.push_back(trees[0]); 
  if(not runParams.isHeatedChainsUseSame())
    for(auto iter =  trees.begin() + 1  ; iter < trees.end(); ++iter)
      treesToInitialize.push_back(*iter); 

  // choose how to initialize the topology
  for(auto &tralnPtr : treesToInitialize)
    {
      auto tr = tralnPtr->getTr();
      bool hasBranchLength = false; 
      if(treesConsumed < numTreesAvailable) // use a user defined starting tree 
	{
	  hasBranchLength = readTreeWithOrWithoutBL(tr, treeFH);
	  ++treesConsumed; 
	}
      else
	{	      
	  if(runParams.isUseParsimonyStarting())
	    TreeRandomizer::createParsimonyTree(*tralnPtr, treeRandomness); 
	  else
	    TreeRandomizer::randomizeTree(*tralnPtr, treeRandomness); 
	}	  
	  
      // correct branches 
      for(auto &b : tralnPtr->extractBranches())
	{
	  b.setLength( hasBranchLength ? 
		       ( - exp(b.getLength() / tr->fracchange)) 
		       : TreeAln::initBL ); 
	  b.applyToTree(*tralnPtr); 
	}
    }

  // cout << "ref tree now is "<< *(trees[0]) << endl; 
  
  // propagate the tree to the coupled chains, if necessar y
  if(runParams.isHeatedChainsUseSame())
    {
      TreeAln &ref = *( trees[0]); 
      bool isFirst = true; 
      for(shared_ptr<TreeAln> &treePtr : trees)
	{
	  if(isFirst )
	    isFirst = false; 
	  else 
	    (*treePtr) = ref; 
	}
    }
}


void SampleMaster::initializeRuns( )
{
  FILE *treeFH = NULL; 
  nat numTreesAvailable = countNumberOfTreesQuick(cl.getTreeFile().c_str()); 

  // initialize one tree 
  vector<shared_ptr<TreeAln> > trees; 
  trees.push_back(make_shared<TreeAln>());
  trees[0]->initializeFromByteFile(cl.getAlnFileName()); 
  trees[0]->enableParsimony();

  vector<unique_ptr<AbstractProposal> > proposals; 
  vector<shared_ptr<RandomVariable> > variables; 

  auto restorer = make_shared<LnlRestorer>(*(trees[0]));
  auto eval = make_shared<LikelihoodEvaluator>(restorer); 
  
  initWithConfigFile(cl.getConfigFileName(), trees[0], proposals, variables, eval);
  assert(runParams.getTuneFreq() > 0); 
  
#ifdef DEBUG_LNL_VERIFY
  auto dT = make_shared<TreeAln>();
  dT->initializeFromByteFile(cl.getAlnFileName()); 
  dT->enableParsimony(); 
  eval->setDebugTraln(dT);
#endif

  // ORDER: must be after  initWithConfigFile
  tout << endl << "Will run " << runParams.getNumRunConv() << " runs in total with "<<  runParams.getNumCoupledChains() << " coupled chains" << endl << endl;   
  for(int i = trees.size(); i < runParams.getNumCoupledChains(); ++i)
    {
      auto traln =  make_shared<TreeAln>();
      trees.push_back(traln); 
      trees[i]->initializeFromByteFile(cl.getAlnFileName()); 
      trees[i]->enableParsimony();
    }

  if( numTreesAvailable > 0 )
    treeFH = myfopen(cl.getTreeFile().c_str(), "r"); 

  vector<randCtr_t> runSeeds;
  vector<randCtr_t> treeSeeds; 
  for(int i = 0; i < runParams.getNumRunConv();++i)
    {
      runSeeds.push_back(masterRand.generateSeed()); 
      treeSeeds.push_back(masterRand.generateSeed()); 
    }

  nat treesConsumed = 0; 
  for(int i = 0; i < runParams.getNumRunConv() ; ++i)
    {    
      initTrees(trees, treeSeeds[i], treesConsumed, numTreesAvailable, treeFH); 

      if( i % pl.getRunsParallel() == pl.getMyRunBatch() )
	{
	  vector<Chain> chains; 

	  for(int j = 0; j < runParams.getNumCoupledChains(); ++j)
	    {
	      randCtr_t c; 
	      chains.push_back(Chain( c , trees[j], proposals, eval ) ); 
	      auto &chain = chains[j]; 		
	      chain.setRunId(i); 
	      chain.setTuneFreuqency(runParams.getTuneFreq()); 
	      chain.setHeatIncrement(j); 
	      chain.setDeltaT(runParams.getHeatFactor()); 
	    }
	  
	  CoupledChains run(runSeeds[i], i, cl.getWorkdir(), runParams.getNumCoupledChains(), chains) ; 
	  run.setTemperature(runParams.getHeatFactor());
	  run.setTuneHeat(runParams.getTuneHeat()); 
	  run.setPrintFreq(runParams.getPrintFreq()); 
	  run.setSwapInterval(runParams.getSwapInterval()); 
	  run.setSamplingFreq(runParams.getSamplingFreq()); 
	  run.setRunName(runParams.getRunId()); 
	  run.seedChains(); 
	  runs.push_back(run); 	  
	  
	  cout << runs[runs.size()-1].getChains()[0].getTraln() << endl; 
	}
      cout << "================================================================"<< endl;
    }


  // determine and print initial state for each chain  
  tout << "initial state: " << endl; 
  tout << "================================================================" << endl; 
  for(auto &run: runs)
    {
      bool isFist = true; 
      for(auto &chain :  run.getChains())
	{
	  if(isFist)
	    {
	      cout << chain.getTraln() << endl; 
	      isFist = false ; 
	    }

	  chain.resume(); 
	  tout << chain; 
	  tout << "\tRNG(" << chain.getChainRand() << ")"<< endl; 	  
	  // tout << chain.getTraln() << endl; 
	}
      tout << "================================================================" << endl; 
    }
  
  if(numTreesAvailable > 0)
     fclose(treeFH); 

  tout << endl; 
}


//STAY 
bool SampleMaster::convergenceDiagnostic()
{
  if(runParams.getNumRunConv() > 1)    
    { 
      vector<string> fns; 
      for(int i = 0; i < runParams.getNumRunConv(); ++i)
	{
	  stringstream ss; 
	  ss <<  PROGRAM_NAME << "_topologies." << runParams.getRunId() << "." << i; 
	  fns.push_back(ss.str());
	}
     
      AvgSplitFreqAssessor asdsf(fns);

      int end = asdsf.getEnd();
      
      int treesInBatch = runParams.getDiagFreq() / runParams.getSamplingFreq(); 

      end /= treesInBatch; 
      end *= treesInBatch;       

      if(end > 0)
	{	  
	  asdsf.setEnd(end);
	  if( runParams.getBurninGen() > 0 )
	    {
	      assert(runParams.getBurninProportion() == 0.); 

	      int treesToDiscard =  runParams.getBurninGen() / runParams.getSamplingFreq(); 

	      if(end < treesToDiscard + 2 )
		return false; 
	      else 
		asdsf.setStart(treesToDiscard);  
	    }
	  else 
	    {
	      assert(runParams.getBurninGen() == 0); 
	      int start = (int)((double)end * runParams.getBurninProportion()  ); 
	      asdsf.setStart(start);
	    } 

	  asdsf.extractBips();
	  double asdsfVal = asdsf.computeAsdsf(runParams.getAsdsfIgnoreFreq());

#if HAVE_PLL == 0      
	  if(processID == 0)
#endif 
	    tout  << "ASDSF for trees " << asdsf.getStart() << "-" << asdsf.getEnd() << ": " <<setprecision(2) << asdsfVal * 100 << "%" << endl; 

	  return asdsfVal < runParams.getAsdsfConvergence(); 

	}
      else 
	return false; 
      
    }
  else 
    {
      return runs[0].getChains()[0].getGeneration() > runParams.getNumGen();
    }
}


// static void initWithStartingTree(FILE *fh, vector<shared_ptr<TreeAln> > &tralns)
// {
//   // fetch a tree 
//   auto traln = tralns[0]; 
//   tree *tr = traln->getTr();


//   int count = 0;       


//   for(auto &b : extractBranches )

//   // if(hasBranchLength)


//   if(hasBranchLength)
//     traverseInitCorrect(tr->start->back, &count, traln ) ;  
//   else      
//     traverseInitFixedBL(tr->start->back, &count, traln, TreeAln::initBL ); 

//   assert(count == 2 * tr->mxtips  -3);       

//   for(nat i = 1 ; i < tralns.size(); ++i)
//     *(tralns[i]) = *traln; 
// }


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


void SampleMaster::initWithConfigFile(string configFileName, shared_ptr<TreeAln> traln, 
				      vector<unique_ptr<AbstractProposal> > &proposalResult, 
				      vector<shared_ptr<RandomVariable> > &variableResult, shared_ptr<LikelihoodEvaluator> eval)
{
  ConfigReader reader; 
  ifstream fh(configFileName); 
  NxsToken token(fh); 

  paramBlock.setTree(traln); 
  reader.Add(&paramBlock); 

  BlockPrior priorBlock(traln->getNumberOfPartitions());
  reader.Add(&priorBlock);

  reader.Add(&runParams);

  BlockProposalConfig proposalConfig; 
  reader.Add(&proposalConfig);   

  reader.Execute(token);

  RunFactory r; 
  r.configureRuns(proposalConfig, priorBlock , paramBlock, *traln, proposalResult, eval);
  variableResult = r.getRandomVariables(); 
}


void SampleMaster::validateRunParams()
{
  assert(runParams.getNumCoupledChains() > 0); 
  // TODO 
}


// a developmental mode to integrate over branch lengths
// #define _GO_TO_INTEGRATION_MODE


void SampleMaster::run()
{
  bool hasConverged = false;   
  while(not hasConverged)   
    {      
      for(auto &run : runs)
	run.executePart(runParams.getDiagFreq());

      hasConverged = convergenceDiagnostic(); 

#ifdef ENABLE_PRSF      
      printPRSF(run_id);
#endif
    }

#ifdef _GO_TO_INTEGRATION_MODE
  // go into integration mode, if we want to do so 
  branchLengthsIntegration();
#endif
}
 
 
void SampleMaster::finalizeRuns()
{
  for(auto &run : runs)
    {
      auto &chains = run.getChains(); 
      for(auto &chain : chains)
	{
	  if( chain.getChainHeat() == 1.f)
	    chain.finalizeOutputFiles(run.getTopoFile());
	  tout << "best state was: " << chain.getBestState( )<< endl; 
	}
    }

  tout << endl << "Converged/stopped after " << runs[0].getChains()[0].getGeneration() << " generations" << endl;   
  tout << endl << "Total execution time: " 
       << CLOCK::duration_cast<CLOCK::duration<double> >( CLOCK::system_clock::now() - initTime   ).count() <<  " seconds" << endl; 
}


#include "IntegrationModuleImpl.hpp"
