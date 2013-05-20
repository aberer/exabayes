#include <sstream>
#include <fstream>
#include <memory>


#include "BlockRunParameters.hpp"
#include "ConfigReader.hpp"
#include "output.h"
#include "eval.h"
#include "SampleMaster.hpp"
#include "Chain.hpp"
#include "TreeRandomizer.hpp"
#include "treeRead.h"
#include "tune.h"
#include "ExtendedTBR.hpp"
#include "ExtendedSPR.hpp"
#include "ParsimonySPR.hpp"
#include "StatNNI.hpp"
#include "RadiusMlSPR.hpp"
#include "Category.hpp"
#include "NodeSlider.hpp"
#include "BranchLengthMultiplier.hpp"
#include "BranchCollapser.hpp"
#include "AminoModelJump.hpp"


// TODO =( 
#include "GlobalVariables.hpp"

#include "LnlRestorer.hpp"
#include "AvgSplitFreqAssessor.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"
#include "TreeLengthMultiplier.hpp"

#include "PartitionProposal.hpp"

extern double masterTime; 

static int countNumberOfTreesQuick(const char *fn ); 
static void initWithStartingTree(FILE *fh, vector<TreeAln*> &tralns); 
static void initTreeWithOneRandom(int seed, vector<TreeAln*> &tralns); 

SampleMaster::SampleMaster(const CommandLine &cl , const ParallelSetup &pl) 
  : esprStopProp(0.5)
  , parsimonyWarp(0.1)
  , guidedRadius(5)
  , pl(pl)
  , initTime(gettime())
{
  FILE *treeFH = NULL; 
  int numTreesAvailable = countNumberOfTreesQuick(cl.getTreeFile().c_str()); 

  PriorBelief prior;
  vector<Category> proposals; 
  vector<double> proposalWeights; 

  vector<TreeAln*> trees; 
  initTrees(trees, cl );

  initWithConfigFile(cl.getConfigFileName(), prior, proposalWeights, *(trees[0]));
  assert(runParams.getTuneFreq() > 0); 
  assert(esprStopProp > 0) ;
  setupProposals(proposals, proposalWeights, prior);

  if( numTreesAvailable > 0 )
    treeFH = myfopen(cl.getTreeFile().c_str(), "r"); 
  
  Randomness masterRand(cl.getSeed());   
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
      if( i < numTreesAvailable)
	initWithStartingTree(treeFH, trees); 
      else 
	initTreeWithOneRandom(treeSeeds[i], trees);

      if( i %  pl.getRunsParallel() == pl.getMyRunBatch() )
	runs.push_back(CoupledChains(runSeeds[i], runParams.getNumCoupledChains(), trees, i, runParams.getPrintFreq(), runParams.getSwapInterval(), runParams.getSamplingFreq(), runParams.getHeatFactor(), cl.getRunid(), cl.getWorkdir(), prior, proposals, runParams.getTuneFreq())); 
    }
  
	  if(runParams.getTuneHeat())
    for(CoupledChains& r : runs)
      r.enableHeatTuning(runParams.getTuneFreq()); 

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
      return runs[0].getChain(0)->getGeneration() > runParams.getNumGen();
    }
}


// LEGACY
static void traverseInitFixedBL(nodeptr p, int *count, TreeAln *traln,  double z )
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


void SampleMaster::initWithConfigFile(string configFileName, PriorBelief &prior, vector<double> &proposalWeights, const TreeAln &traln )
{
  ConfigReader reader; 
  ifstream fh(configFileName); 
  NxsToken token(fh); 

  BlockParams paramBlock(traln.getNumberOfPartitions()); 
  reader.Add(&paramBlock); 

  BlockPrior priorBlock(traln.getNumberOfPartitions());
  reader.Add(&priorBlock);

  reader.Add(&runParams);

  reader.Execute(token);

  // auto parameters =  paramBlock.getParameters() ;
  // cout << "parsed parameters: " ; 
  // for(auto v : parameters)
  //   cout << v << ","; 
  // cout << endl; 
  // assert(0); 

  
  // ConfigHandler conf(configFileName); 

  assert(NOT_IMPLEMENTED); 

}


void SampleMaster::validateRunParams()
{
  assert(runParams.getNumCoupledChains() > 0); 
  // TODO 
}


void SampleMaster::initTrees(vector<TreeAln*> &trees, const CommandLine &cl )
{  
  tout << endl << "Will run " << runParams.getNumRunConv() << " runs in total with "<<  runParams.getNumCoupledChains() << " coupled chains" << endl << endl; 

#ifdef DEBUG_LNL_VERIFY
  globals.debugTree = new TreeAln();   
  globals.debugTree->initializeFromByteFile(cl.getAlnFileName()); 
  globals.debugTree->enableParsimony();
#endif

  for(int i = 0; i < runParams.getNumCoupledChains(); ++i)
    {
      TreeAln *traln = new TreeAln();
      traln->initializeFromByteFile(cl.getAlnFileName()); 
#if HAVE_PLL != 0 
      traln->enableParsimony();
#endif
      trees.push_back(traln); 
    }
  
  // only use one restorer for all chains 
  auto restorer = shared_ptr<LnlRestorer> (new LnlRestorer(*(trees[0])));
  for(auto tree : trees)
    tree->setRestorer(restorer);
}


void SampleMaster::run()
{
  bool hasConverged = false;   
  while(not hasConverged)   
    {      
      for(nat i = 0; i < runs.size(); ++i)
	{
	  auto run = runs[i]; 
	  run.executePart(runParams.getDiagFreq());
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
		chain->finalizeOutputFiles(run.getTopoFile());
	    }
	}

      tout << endl << "Converged after " << runs[0].getChain(0)->getGeneration() << " generations" << endl; 
      tout << endl << "Total execution time: " << setprecision(2) <<  gettime() - initTime <<  " seconds" << endl; 
    }  
}


/** 
    @brief sets up the proposals depending on the prior configuration 
    
 */ 
void SampleMaster::setupProposals(vector<Category> &proposalCategories, vector<double> proposalWeights, const PriorBelief &prior)
{
  vector<AbstractProposal*> prop; 
  
  vector<aaMatrix_t> someMatrices; // TODO, probably extract from prior belief 

  // initialize proposals 
  for(int i = 0; i < NUM_PROPOSALS ; ++i)
    { 
      double weight = proposalWeights[proposal_type(i)]; 
      if( weight != 0)
	{
	  AbstractProposal *proposal = NULL; 
	  
	  switch(proposal_type(i))
	    {	      
	    case BRANCH_LENGTHS_MULTIPLIER:	      
	      proposal = new BranchLengthMultiplier(weight, INIT_BL_MULT) ; 
	      break; 
	    case NODE_SLIDER:
	      proposal = new NodeSlider( weight, INIT_NODE_SLIDER_MULT); 
	      break; 
	    case REVMAT_SLIDER: 
	      proposal = new PartitionProposal<SlidingProposal, RevMatParameter>(weight, INIT_RATE_SLID_WIN, "revMatSlider"); 
	      break; 
	    case FREQUENCY_SLIDER:
	      proposal = new PartitionProposal<SlidingProposal, FrequencyParameter>( weight, INIT_FREQ_SLID_WIN, "freqSlider"); 
	      break; 		  
	    case TL_MULT:
	      proposal = new TreeLengthMultiplier( weight, INIT_TL_MULTI); 
	      break; 
	    case E_TBR: 
	      proposal = new ExtendedTBR( weight, esprStopProp, INIT_ESPR_MULT); 
	      break; 
	    case E_SPR: 
	      proposal = new ExtendedSPR( weight, esprStopProp, INIT_ESPR_MULT); 
	      break; 
	    case PARSIMONY_SPR:	
	      proposal = new ParsimonySPR( weight, parsimonyWarp, INIT_ESPR_MULT); 
	      break; 
	    case ST_NNI: 
	      proposal = new StatNNI( weight, INIT_NNI_MULT); 
	      break; 
	    case RATE_HET_MULTI: 
	      proposal = new PartitionProposal<MultiplierProposal,RateHetParameter>( weight, INIT_GAMMA_MULTI, "rateHetMulti"); 
	      break; 
	    case RATE_HET_SLIDER: 
	      proposal = new PartitionProposal<SlidingProposal,RateHetParameter>( weight, INIT_GAMMA_SLID_WIN, "rateHetSlider"); 
	      break; 
	    case FREQUENCY_DIRICHLET: 
	      proposal = new PartitionProposal<DirichletProposal,FrequencyParameter>( weight, INIT_DIRICHLET_ALPHA, "freqDirich"); 
	      break; 
	    case REVMAT_DIRICHLET: 
	      proposal = new PartitionProposal<DirichletProposal,RevMatParameter>(weight, INIT_DIRICHLET_ALPHA, "revMatDirich"); 	      
	      break; 
	    case GUIDED_SPR:
	      proposal = new RadiusMlSPR( weight, guidedRadius ); 
	      break; 
	    case BRANCH_COLLAPSER:
	      proposal = new BranchCollapser(weight); 
	      break; 
	    case AMINO_MODEL_JUMP: 
	      proposal = new AminoModelJump(weight, someMatrices);
	      break; 
	    default : 
	      assert(0); 
	    }

	  if(not prior.categoryIsFixed(proposal->getCategory()))
	    prop.push_back(proposal);
	  else 
	    tout << "Notice: Relative weight " << weight << " was specified for proposal "   << proposal->getName()  << ". Since the prior for this category was set to fixed, this proposal will be ignored." << endl; 	      
	}
    }


  // get total sum 
  double sum = 0; 
  for(auto p : prop)
    sum += p->getRelativeProbability();

  // create categories 
  vector<string> allNames = {"Topology", "BranchLengths", "Frequencies", "RevMatrix", "RateHet" }; 
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    {
      // fish out the correct proposals 
      vector<AbstractProposal*> pr; 
      double catSum = 0; 
      for(auto p : prop)
	{	  
	  if(p->getCategory() == i)
	    {
	      pr.push_back(p); 
	      catSum += p->getRelativeProbability();
	    }
	}

      if(pr.size() > 0)
	proposalCategories.push_back( Category(allNames[i], category_t(i), catSum / sum, pr )); 
    }    


  if ( isOutputProcess() )
    {
      // print some info 
      tout << "using the following moves: " << endl; 
      for(nat i = 0; i < proposalCategories.size(); ++i)
	{
	  // TODO also to info file 
      
	  tout << proposalCategories[i].getName() << " " << fixed << setprecision(2) << proposalCategories[i].getCatFreq() * 100 << "%\t" ; 
	  auto cat = proposalCategories[i]; 
	  auto p1 =  cat.getProposals(); 
	  for(auto p : p1)
	    tout << "\t" << p->getName() << "(" << fixed << setprecision(2) <<  p->getRelativeProbability() * 100 << "%)" ;
	  tout << endl; 
	}
    }
}
