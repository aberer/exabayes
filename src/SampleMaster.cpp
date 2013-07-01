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
static void initWithStartingTree(FILE *fh, vector<TreeAlnPtr > tralns); 
static void initTreeWithOneRandom(int seed,  vector<TreeAlnPtr > &tralns); 

SampleMaster::SampleMaster(const ParallelSetup &pl, const CommandLine& _cl ) 
  : pl(pl)
  , initTime(CLOCK::system_clock::now())
  , masterRand(cl.getSeed())
  , cl(_cl)
{
  initializeRuns(); 
}



void SampleMaster::initializeRuns( )
{
  FILE *treeFH = NULL; 
  int numTreesAvailable = countNumberOfTreesQuick(cl.getTreeFile().c_str()); 

  // initialize one tree 
  vector<TreeAlnPtr> trees; 
  trees.push_back(TreeAlnPtr(new TreeAln()));
  trees[0]->initializeFromByteFile(cl.getAlnFileName()); 
  trees[0]->enableParsimony();

  vector<ProposalPtr> proposals; 
  vector<RandomVariablePtr> variables; 

  LnlRestorerPtr restorer(new LnlRestorer(*(trees[0])));
  LikelihoodEvaluatorPtr eval(new LikelihoodEvaluator(restorer));

  initWithConfigFile(cl.getConfigFileName(), trees[0], proposals, variables, eval);
  assert(runParams.getTuneFreq() > 0); 

  // ORDER: must be after  initWithConfigFile
  initTrees(trees);

  if( numTreesAvailable > 0 )
    treeFH = myfopen(cl.getTreeFile().c_str(), "r"); 

  vector<int> runSeeds; 
  vector<int> treeSeeds; 
  for(int i = 0; i < runParams.getNumRunConv();++i)
    {
      randCtr_t r = masterRand.generateSeed();
      runSeeds.push_back(r.v[0]); 
      treeSeeds.push_back(r.v[1]); 
    }

  for(int i = 0; i < runParams.getNumRunConv() ; ++i)
    {      
      assert(i < 1); 
      if( i < numTreesAvailable)
	initWithStartingTree(treeFH, trees); 
      else 
	initTreeWithOneRandom(treeSeeds[i], trees);

      if( i %  pl.getRunsParallel() == pl.getMyRunBatch() )
	runs.push_back(CoupledChainsPtr(new CoupledChains(runSeeds[i], i, runParams, trees, cl.getWorkdir(), proposals, variables, eval))); 
    }

  if(runParams.getTuneHeat())
    {
      for(CoupledChainsPtr& r : runs)
	r->enableHeatTuning(runParams.getTuneFreq()); 
    }

  if(numTreesAvailable > 0)
     fclose(treeFH); 

  cout << "we have " << trees[0]->getNumberOfPartitions() << " partitions"  << endl; 

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
      return runs[0]->getChains()[0]->getGeneration() > runParams.getNumGen();
    }
}


// LEGACY
static void traverseInitFixedBL(nodeptr p, int *count, shared_ptr<TreeAln> traln,  double z )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;
  
  for( i = 0; i < traln->getNumBranches(); i++)
      traln->setBranchLengthBounded(z, i, p); 
  
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


static void initTreeWithOneRandom(int seed,  vector<TreeAlnPtr> &tralns)
{  
  TreeRandomizer trRandomizer(seed);
  for(auto &traln : tralns)
    {
      trRandomizer.randomizeTree(*traln); 
      int count = 0; 
      traverseInitFixedBL( traln->getTr()->start->back, &count, traln, TreeAln::initBL);
      assert(count  == 2 * traln->getTr()->mxtips - 3);
    }
}



static void initWithStartingTree(FILE *fh, vector<TreeAlnPtr> tralns)
{
  // fetch a tree 
  auto traln = tralns[0]; 
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


void SampleMaster::initWithConfigFile(string configFileName, const TreeAlnPtr &traln, vector<ProposalPtr> &proposalResult, 
				      vector<RandomVariablePtr> &variableResult, LikelihoodEvaluatorPtr &eval)
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



void SampleMaster::initTrees(vector<TreeAlnPtr> &trees )
{  
  tout << endl << "Will run " << runParams.getNumRunConv() << " runs in total with "<<  runParams.getNumCoupledChains() << " coupled chains" << endl << endl; 
#ifdef DEBUG_LNL_VERIFY
  globals.debugTree = new TreeAln();   
  globals.debugTree->initializeFromByteFile(cl.getAlnFileName()); 
  globals.debugTree->enableParsimony();
#endif

  for(int i = trees.size(); i < runParams.getNumCoupledChains(); ++i)
    {
      auto traln =  TreeAlnPtr(new TreeAln);
      traln->initializeFromByteFile(cl.getAlnFileName()); 
      traln->enableParsimony();
      trees.push_back(traln); 
    }

}


// a developmental mode to integrate over branch lengths
// #define _GO_TO_INTEGRATION_MODE


void SampleMaster::run()
{
  bool hasConverged = false;   
  while(not hasConverged)   
    {      
      for(auto &run : runs)
	run->executePart(runParams.getDiagFreq());

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
      auto chains = run->getChains(); 
      for(auto &chain : chains)
	{
	  if( chain->getChainHeat() == 1.f)
	    chain->finalizeOutputFiles(run->getTopoFile());
	  tout << "best state was: " << chain->getBestState( )<< endl; 
	}
    }

  tout << endl << "Converged/stopped after " << runs[0]->getChains()[0]->getGeneration() << " generations" << endl;   
  tout << endl << "Total execution time: " 
       << CLOCK::duration_cast<CLOCK::duration<double> >( CLOCK::system_clock::now() - initTime   ).count() <<  " seconds" << endl; 
}


#define STEPS_FOR_LNL 1000
#define INTEGRATION_GENERATIONS 10000
#define NR_STEPS 30

#include <sstream>
#include "proposals/BranchIntegrator.hpp"
#include "ProposalRegistry.hpp"

void SampleMaster::branchLengthsIntegration()  
{
  assert(runs.size() == 1 );   
  auto &run = runs[0];   
  auto chains = run->getChains(); 
  assert(chains.size() == 1); 
  auto chain = chains[0]; 
  
  stringstream ss; 
  ss << cl.getRunid() << ".tree.tre" ; 

  ofstream tFile( ss.str());   
  TreePrinter tp(true, true, false);   
  TreePrinter tp2(true, true, true); 
  auto traln = chain->getTraln(); 
  tFile << tp.printTree(traln) << endl; 
  tFile << tp2.printTree(traln) << endl; 
  tFile.close(); 
  
  auto eval = chain->getEvaluator();

  RandomVariablePtr r(new RandomVariable(Category::BRANCH_LENGTHS, 0));
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    r->addPartition(i);

  vector<RandomVariablePtr> vars =  {r} ; 

  double lambda   =  10 ; 


  vars[0]->setPrior(PriorPtr (new ExponentialPrior(lambda)));

  auto p = ProposalPtr(new BranchIntegrator(ProposalRegistry::initBranchLengthMultiplier));   
  vector<ProposalPtr>  proposals;   
  proposals.push_back( std::move(p) ); 
  proposals[0]->addPrimVar(vars[0]); 

  Chain integrationChain(masterRand.generateSeed(), 
  			 0, 
  			 0, 
  			 TreeAlnPtr(&(traln)), 
  			 proposals, 
  			 100,
  			 vars, 
  			 eval ); 

  auto branches =  traln.extractBranches();
  vector<ProposalPtr>& ps = integrationChain.getProposals(); 
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0].get()); 

  int ctr = 0; 
  for(auto branch : branches)
    {      
      double minHere = 1000; 
      double maxHere = 0; 

      branch.updateLength(traln); 
      
      integrator->setToPropose(branch); 
      
      stringstream ss; 
      ss << "samples." << cl.getRunid()<< "." << branch.getPrimNode() << "-" << branch.getSecNode()   <<  ".tab" ;
      ofstream thisOut (ss.str()); 
      
      // run the chain to integrate 
      cout << "integrating branch " << branch << endl; 
      for(int i = 0; i < INTEGRATION_GENERATIONS; ++i) 
	{	  
	  integrationChain.step();
	  auto p = branch.findNodePtr(traln); 	  
	  branch.updateLength(traln); 

	  double length = branch.getInterpretedLength(traln); 
	  thisOut << length << endl; 
	  
	  if(length < minHere)
	    minHere = length; 
	  if(maxHere < length )
	    maxHere = length; 
	}      

      thisOut.close();
      
      // get branch lengths 
      ss.str(std::string()); 
      ss << "lnl." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str());

      cout << "evaluating branch lengths for " << branch << endl; 
      
      if(maxHere != minHere)
	{
	  for(double i = minHere; i < maxHere+0.00000001 ; i+= (maxHere-minHere)/ STEPS_FOR_LNL)
	    {
	      double tmp = branch.getInternalLength(traln,i); 
	      traln.setBranchLengthBounded(tmp, 0, branch.findNodePtr(traln)); 
	      eval->evaluate(traln, branch, false);
	      double lnl = traln.getTr()->likelihood; 
	      
	      thisOut << i << "\t" << setprecision(std::numeric_limits<double>::digits10) << lnl << endl; 
	    }
	}
      else
	thisOut << minHere << "\t" << "NA" << endl; 
      thisOut.close(); 

      Branch tmpBranch = branch; 
      cout << "optimizing the branch using nr" << endl; 
      ss.str(std::string());
      ss << "nr-length." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str()); 

      double result = 0;  
      double curVal = branch.getInternalLength(traln,0.1); 
      double secDerivative = 0; 
      double firstDerivative = 0; 

      double prevVal = curVal; 
      for(int i = 0; i < NR_STEPS; ++i )
	{
	  makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), 
			  branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
			  &curVal, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
	  tmpBranch.setLength(result);
	  thisOut << prevVal <<  "\t" << firstDerivative << "\t" << secDerivative << endl; 	
	  prevVal = tmpBranch.getInterpretedLength(traln); 
	  curVal = result; 
	} 

      double something = tmpBranch.getInternalLength(traln, prevVal); 

      makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
      
      thisOut << prevVal << "\t" << firstDerivative << "\t" << secDerivative << endl; 

      thisOut.close(); 
      
      // reset 
      double tmp = branch.getLength();      
      traln.setBranchLengthBounded(tmp , 0,branch.findNodePtr(traln) );
    }

  cout << "finished!" << endl; 
}
