#include <sstream>
#include <fstream>
#include <memory>

#include "AdHocIntegrator.hpp"
#include "ProposalSet.hpp"

#include "common.h"

#include "config/BlockProposalConfig.hpp"
#include "config/BlockRunParameters.hpp"
#include "config/ConfigReader.hpp"

#include "priors/FixedPrior.hpp"
#include "parameters/TopologyParameter.hpp"

#include "config/MemoryMode.hpp"
#include "SampleMaster.hpp"
#include "Chain.hpp"
#include "TreeRandomizer.hpp"
#include "treeRead.h"
#include "tune.h"
#include "RunFactory.hpp"	

#include "GlobalVariables.hpp"


#include "BoundsChecker.hpp"
#include "ArrayRestorer.hpp"
#include "AvgSplitFreqAssessor.hpp"

#include "proposers/AbstractProposer.hpp"
#include "proposers/SlidingProposal.hpp"
#include "proposers/MultiplierProposal.hpp"

#include "PlainLikelihoodEvaluator.hpp"
#include "RestoringLnlEvaluator.hpp"

#include "NewRestoringEvaluator.hpp"

#include "TreeIntegrator.hpp"

// TODO 
// * check restart 

// a developmental mode to integrate over branch lengths

// #define _GO_TO_INTEGRATION_MODE


// #ifndef _EXPERIMENTAL_INTEGRATION_MODE
// #undef _GO_TO_INTEGRATION_MODE
// #endif

void genericExit(int code); 


SampleMaster::SampleMaster(const ParallelSetup &pl, const CommandLine& _cl ) 
  : pl(pl)
  , initTime(CLOCK::system_clock::now())
  , cl(_cl)
  , lastPrintTime(CLOCK::system_clock::now())
  , masterRand(cl.getSeed())
{
}


void SampleMaster::initializeTree(TreeAln &traln, std::string startingTree, Randomness &treeRandomness, const std::vector<AbstractParameter*> &params)
{  
  auto tr = traln.getTr();
  bool hasBranchLength = false; 
  if(startingTree.compare("") != 0 ) // meh 
    hasBranchLength = readTreeWithOrWithoutBL(tr, startingTree);
  else
    {	      
      if(runParams.isUseParsimonyStarting())
	TreeRandomizer::createParsimonyTree(traln, treeRandomness); 
      else
	TreeRandomizer::randomizeTree(traln, treeRandomness); 
    }

  for(auto &b : traln.extractBranches())
    {
      for(auto &param : params)
	{
	  double initLen =  TreeAln::initBL ; 

	  if(hasBranchLength)
	    initLen = traln.getBranch(b,param).getLength() ; 

	  auto bl = b.toBlDummy(); 
	  bl.setConvertedInternalLength(traln, param, initLen); 

	  if(not BoundsChecker::checkBranch(bl))
	    {
	      tout << "WARNING: initial branch length already out of bounds. Relative representation: "
		   << bl.getLength() << "\treal length=" << bl.getInterpretedLength(traln, param) << std::endl; 
	      BoundsChecker::correctBranch(bl); 	      
	    }
	  traln.setBranch(bl, param); 
	}
    }
}



void SampleMaster::initTrees(vector<shared_ptr<TreeAln> > &trees, randCtr_t seed, nat &treesConsumed, std::vector<std::string> startingTreeStrings, const std::vector<AbstractParameter*> &params)
{  
  auto treeRandomness = Randomness(seed); 

  vector<shared_ptr<TreeAln> > treesToInitialize; 
  treesToInitialize.push_back(trees[0]); 
  if(not runParams.isHeatedChainsUseSame())
    for(auto iter =  trees.begin() + 1  ; iter < trees.end(); ++iter)
      treesToInitialize.push_back(*iter); 

  // choose how to initialize the topology
  for(auto &tralnPtr : treesToInitialize)
    {
      auto stTr = 
	treesConsumed < startingTreeStrings.size() ? startingTreeStrings.at(treesConsumed) : std::string{""}; 
      ++treesConsumed; 		// TODO re-engineer!
      initializeTree(*tralnPtr, stTr, treeRandomness, params); 
    }

  // propagate the tree to the coupled chains, if necessary
  if(runParams.isHeatedChainsUseSame())
    {
      TreeAln &ref = *( trees[0]); 
      bool isFirst = true; 
      for(shared_ptr<TreeAln> &treePtr : trees)
	{
	  if(isFirst )
	    isFirst = false; 
	  else 
	    treePtr->copyModel(ref); 
	}
    }
}


void SampleMaster::printAlignmentInfo(const TreeAln &traln)
{
  tout << std::endl << "The (binary) alignment file you provided, consists of the following\n" 
       << "partitions:" 
       << std::endl;
  
  for(nat i = 0 ;i < traln.getNumberOfPartitions() ;++i)
    {
      auto partition = traln.getPartition(i);       
      nat length = partition->upper - partition->lower; 
      
      tout << std::endl; 

      tout << "number:\t\t" << i << std::endl;     
      tout << "name:\t\t" << partition->partitionName << std::endl; 
      tout << "#patterns:\t" << length << std::endl;       

      switch(partition->dataType)
	{
	case DNA_DATA: 
	  tout << "type:\t\tDNA" << std::endl; 
	  break; 
	case AA_DATA: 
	  tout << "type:\t\tAA" << std::endl; 
	  break; 
	default : 
	  assert(0); 
	}      
    }
  tout << "\nParameter names will reference the above partition numbers."
       << std::endl
       << "================================================================"
       << std::endl  ; 
}


void SampleMaster::informPrint()
{
#if HAVE_PLL == 0
  if(pl.getRunsParallel()  > 1 )
    {
      tout << "Will execute "<< pl.getRunsParallel() << " runs in parallel." << std::endl;       
      if(nat (runParams.getNumRunConv()) < pl.getRunsParallel())
	{
	  tout 
	    << "Error: in the configuration file you specified to run " <<  runParams.getNumRunConv() << " independent\n" 
	    << "runs. Your command line indicates, that you want to run " << pl.getRunsParallel() << " of them\n"
	    << "in parallel. To shield you from surprises, " << PROGRAM_NAME << " will conservatively abort."<< std::endl; 
	  ParallelSetup::genericExit(-1); 
	}
    }
  if(pl.getChainsParallel() > 1)
    {
      tout << "Will execute " << pl.getChainsParallel() << " chains in parallel."<< std::endl; 
      if(nat(runParams.getNumCoupledChains()) < pl.getChainsParallel())	
	{
	  tout
	    << "Error: in the configuration file you specified to run " <<  runParams.getNumCoupledChains() << " coupled\n" 
	    << "chains per run. Your command line indicates, that you want to run " << pl.getChainsParallel() << " of them\n"
	    << "in parallel. To shield you from surprises, " << PROGRAM_NAME << " will conservatively abort."<< std::endl; 
	  ParallelSetup::genericExit(-1); 
	}
    }

  tout << std::endl; 
#endif
}



void SampleMaster::initializeFromCheckpoint()
{
  // continue from checkpoint  
  if(cl.getCheckpointId().compare("") != 0)
    {
      auto prevId = cl.getCheckpointId(); 
      
      std::stringstream ss ;       
      ss << cl.getWorkdir() << ( cl.getWorkdir().compare("") == 0  ? "" : "/") 
	 << PROGRAM_NAME << "_checkpoint." << cl.getCheckpointId();       
      auto checkPointFile = ss.str(); 
      std::ifstream chkpnt; 
      Checkpointable::getIfstream(checkPointFile, chkpnt); 
      if( not chkpnt )
	{
	  tout  << "Warning! Could not open / find file " << checkPointFile << std::endl; 
	  tout << "Although extremely unlikely, the previous run may have been killed\n"
	       << "during the writing process of the checkpoint. Will try to recover from backup..." << std::endl; 
	  
	  ss.str(""); 
	  ss << PROGRAM_NAME << "_prevCheckpointBackup." << cl.getCheckpointId();
	  checkPointFile = ss.str();	  
	  Checkpointable::getIfstream(checkPointFile, chkpnt); 
	  
	  if( not chkpnt)
	    {
	      tout << "Could not recover checkpoint file backup either. Probably there is no backup, since not enough generations have completed. Giving up. Please start a new run." << std::endl; 
	      ParallelSetup::genericExit(-1); 
	    }
	  else 
	    tout << "Success! You lost one checkpoint, but we can continue from the backup. " << std::endl; 	  
	}

      readFromCheckpoint(chkpnt); 
      if(pl.isGlobalMaster())
	{
	  for(auto &run : runs)
	    run.regenerateOutputFiles(cl.getWorkdir(), cl.getCheckpointId());
	  nat curGen = runs[0].getChains()[0].getGeneration();
	  diagFile.regenerate(  cl.getWorkdir(), cl.getRunid(),  cl.getCheckpointId(), 
				curGen, runs[0].getSwapInfo().getMatrix().size() * runs.size());
	}
    }
}



// TODO could move this method somewhere else  
std::unique_ptr<LikelihoodEvaluator> 
SampleMaster::createEvaluatorPrototype(const TreeAln &initTree)
{
  auto eval = std::unique_ptr<LikelihoodEvaluator>{}; 
#ifdef _MY_NEW_EVALUATOR
  eval =  std::unique_ptr<NewRestoringEvaluator> { new NewRestoringEvaluator() }; 
#else 
  switch(cl.getMemoryMode())
    {
      
    case MemoryMode::RESTORING: 
      {
	auto restorer = make_shared<ArrayRestorer>(initTree);  
	eval = std::unique_ptr<RestoringLnlEvaluator>( new RestoringLnlEvaluator(restorer)); 
      }
      break; 
    case MemoryMode::PLAIN:
      {
	// NOT IMPLEMENTED YET 
	assert(0); 
	// eval = make_shared<PlainLikelihoodEvaluator>();
      }
      break; 
    default : 
      assert(0);       
    }  
#endif


#ifdef DEBUG_LNL_VERIFY
  auto dT = make_shared<TreeAln>();
  dT->initializeFromByteFile(cl.getAlnFileName(), RunModes::NOTHING ); 
  // even less efficient, but safer 

  dT->enableParsimony(); 
  eval->setDebugTraln(dT);
#endif

  return eval; 
}


void SampleMaster::initializeWithParamInitValues(std::vector<shared_ptr<TreeAln>> &tralns , const std::vector<AbstractParameter*> &params ) const 
{
  for(auto &param : params)
    {
      // no need to do anything for the branch lengths parameter =>
      // that's done via tree initialization
      if(dynamic_cast<BranchLengthsParameter*>(param) != nullptr 
	 ||  dynamic_cast<TopologyParameter*>(param) != nullptr)
	continue; 

      auto&& prior = param->getPrior(); 

      auto content = ParameterContent{}; 
      content = prior->getInitialValue();
	  
      auto &bla = *(tralns[0].get()); 
      param->verifyContent(bla, content); 
      
      for(auto &traln : tralns)
	param->applyParameter(*traln, content); 
      
      auto result = param->extractParameter(*(tralns[0])); 
    }
}



std::vector<std::string> SampleMaster::getStartingTreeStrings()
{
  auto result =  std::vector<std::string>{};

  auto&& ifh = std::ifstream{cl.getTreeFile()}; 
  if(ifh)
    {
      auto aString = std::string{}; 
      while(std::getline(ifh, aString))
	result.push_back(aString); 
    }

  return result; 
}



void SampleMaster::printProposals(const std::vector<unique_ptr<AbstractProposal> > &proposals, const std::vector<ProposalSet> &proposalSets  ) const 
{
  double sum = 0; 
  for(auto &p : proposals)
    sum += p->getRelativeWeight(); 
  for(auto &p : proposalSets)
    sum += p.getRelativeWeight();
  
  tout << "Will employ the following proposal mixture (frequency,type,affected variables): " << endl; 
  for(auto &p : proposals )
    {
      tout << PERC_PRECISION << p->getRelativeWeight() / sum * 100 <<   "%\t" ; 
      p->printShort(tout ) ; 
      tout << endl; 
    }
  if(proposals.size() == 0)
    tout << "None." << std::endl; 

  tout << std::endl; 

  if(runParams.isComponentWiseMH())
    {
      tout << "In addition to that, the following sets below will be executed \n"
	   << "in a sequential manner (for efficiency, see manual for how to\n"
	   << "disable)." << std::endl; 
	for(auto &p : proposalSets )
	  {
	    p.printVerboseAbbreviated(tout, sum);
	    tout << std::endl;       
	  }
    }
}


void SampleMaster::printParameters(const TreeAln &traln, const std::vector<unique_ptr<AbstractParameter> > &params) const 
{
  // merely some printing and we are done  
  tout << std::endl << "Parameters to be integrated (initial values derived from prior): " << std::endl; 
  tout << SOME_FIXED_PRECISION; 
  for(auto &v : params)
    {
      tout << v->getId() << "\t" << v.get()  << std::endl; 
      tout << "\tsub-id:\t" << v->getIdOfMyKind() << std::endl; 
      tout << "\tprior:\t" << v->getPrior() << std::endl; 

      if(dynamic_cast<TopologyParameter*>(v.get())  != nullptr )
	tout << "\tinit value:\tgiven (if tree file available) / random" <<  std::endl; 	
      else if(dynamic_cast<BranchLengthsParameter*>(v.get()) != nullptr)
	tout << "\tinit value:\tgiven (if in tree file) / random" <<  std::endl; 
      else 
	tout << "\tinit value:\t" << v->getPrior()->getInitialValue() << std::endl; 
    }
  tout << "================================================================" << std::endl;
  tout << endl; 
}


void SampleMaster::initializeRuns()
{  
  auto startingTrees = getStartingTreeStrings(); 

  if(cl.getCheckpointId().compare(cl.getRunid()) == 0)
    {
      std::cerr << "You specified >" << cl.getRunid() << "< as runid and intended\n"
		<< "to restart from a previous run with id >" << cl.getCheckpointId() << "<."
		<< "Please specify a new runid for the restart. " << std::endl; 
      ParallelSetup::genericExit(-1); 
    }

  // initialize one tree 
  auto&& initTreePtr = std::unique_ptr<TreeAln>(new TreeAln()); 

  auto runmodes = RunModes::NOTHING; 
  if(cl.isPerPartitionDataDistribution())
    runmodes = runmodes | RunModes::PARTITION_DISTRIBUTION; 

  auto trees =  std::vector<std::shared_ptr<TreeAln> >{}; 
  initTreePtr->initializeFromByteFile(cl.getAlnFileName(), runmodes); 
  initTreePtr->enableParsimony();
  const auto& initTree = *initTreePtr; 

  // START integrator
#ifdef _EXPERIMENTAL_INTEGRATION_MODE
  std::shared_ptr<TreeAln> aTree = std::unique_ptr<TreeAln>(new TreeAln()); 
  aTree->initializeFromByteFile(cl.getAlnFileName(), runmodes); 
  aTree->enableParsimony();

  // let's have another tree for debug
  auto dT = make_shared<TreeAln>();
  dT->initializeFromByteFile(cl.getAlnFileName(), runmodes); 
  dT->enableParsimony(); 

  TreeRandomizer::randomizeTree(*aTree, masterRand); 
  ahInt = new AdHocIntegrator(aTree, dT, masterRand.generateSeed());

  tInt = new TreeIntegrator(aTree, dT, masterRand.generateSeed()); 

#endif
  // END

  auto proposals =  std::vector<unique_ptr<AbstractProposal> >{} ; 
  auto params = std::vector<unique_ptr<AbstractParameter> >{} ; 

  auto evalUptr = createEvaluatorPrototype(initTree); 
  
  auto proposalSets =  std::vector<ProposalSet>{};  
  processConfigFile(cl.getConfigFileName(), &initTree, proposals, params, proposalSets, evalUptr);

  assert(runParams.getTuneFreq() > 0); 

  // ORDER: must be after initWithConfigFile

  for(nat i = 0 ; i < runParams.getNumCoupledChains(); ++i)
    {      
      trees.push_back(make_shared<TreeAln>()); 
      trees[i]->initializeFromByteFile(cl.getAlnFileName(), runmodes); 
      trees[i]->enableParsimony();

#if HAVE_PLL == 0
      if(cl.isPerPartitionDataDistribution())
	{
	  auto tr = trees[i]->getTr(); 
	  tr->manyPartitions = TRUE; 
	}
#endif
    }

  auto runSeeds = vector<randCtr_t>{};
  auto treeSeeds = vector<randCtr_t>{}; 
  for(nat i = 0; i < runParams.getNumRunConv();++i)
    {
      runSeeds.push_back(masterRand.generateSeed()); 
      treeSeeds.push_back(masterRand.generateSeed()); 
    }

  // determine if topology is fixed 
  bool topoIsFixed = false; 
  for(auto &v :params)
    {
      if( dynamic_cast<TopologyParameter*>(v.get())  != nullptr
	  && dynamic_cast<FixedPrior*>(v->getPrior()) != nullptr )
	topoIsFixed = true; 
    }
  
  // gather branch length parameters
  std::vector<AbstractParameter*> blParams;
  for(auto &v : params)
    if(v->getCategory() == Category::BRANCH_LENGTHS)
      blParams.push_back(v.get());

  if(startingTrees.size() > 1 )
    {
      tout << "You provided " << startingTrees.size() << " starting trees" << std::endl;       
      if(topoIsFixed)
	tout << "Since the topology is fixed, only the first one will be used." << std::endl; 
    }

  if(topoIsFixed)
    {
      Randomness treeRandomness(treeSeeds[0]); 
      TreeAln &something = *initTreePtr; 
      initializeTree(something, 
		     startingTrees.size() > 0  ? startingTrees.at(0): std::string(""), // meh 
		     treeRandomness, blParams); 
    }

  nat treesConsumed = 0; 
  for(nat i = 0; i < runParams.getNumRunConv() ; ++i)
    {    
      if(topoIsFixed)
	{
	  for(auto &t : trees)
	    t->copyModel(*initTreePtr); 
	}
      else 	
	initTrees(trees, treeSeeds[i], treesConsumed, startingTrees, blParams); 

      auto paramView = std::vector<AbstractParameter*>{}; 
      for(auto &param : params)
	paramView.push_back(param.get()); 

      initializeWithParamInitValues(trees, paramView );

      auto chains = vector<Chain>{};       
      for(nat j = 0; j < runParams.getNumCoupledChains(); ++j)
	{
	  randCtr_t c; 
	  auto &t = trees.at(j);
	  chains.emplace_back( c , t, proposals, proposalSets, 
			       std::unique_ptr<LikelihoodEvaluator>(evalUptr->clone())  ); 
	  auto &chain = chains[j]; 		
	  chain.setRunId(i); 
	  chain.setTuneFreuqency(runParams.getTuneFreq()); 
	  chain.setHeatIncrement(j); 
	  chain.setDeltaT(runParams.getHeatFactor()); 
	}
	  
      runs.emplace_back(runSeeds[i], i, cl.getWorkdir(), cl.getRunid(), runParams.getNumCoupledChains(), chains); 
      auto &run = *(runs.rbegin()); 
      run.setTemperature(runParams.getHeatFactor());
      run.setTuneHeat(runParams.getTuneHeat()); 
      run.setPrintFreq(runParams.getPrintFreq()); 
      run.setSwapInterval(runParams.getSwapInterval()); 
      run.setSamplingFreq(runParams.getSamplingFreq()); 
      run.seedChains(); 	  
      
      if(pl.isRunLeader() && pl.isMyRun(run.getRunid()))
	run.initializeOutputFiles();
    }

  initializeFromCheckpoint(); 

  if(pl.isGlobalMaster() && not diagFile.isInitialized())
    diagFile.initialize(cl.getWorkdir(), cl.getRunid(), runs);

  // post-pone all printing to the end  
  printAlignmentInfo(initTree);   
  printParameters(initTree, params);
  printProposals(proposals, proposalSets); 
  informPrint();
  printInitializedFiles(); 
  printInitialState(pl); 

  pl.printLoadBalance(initTree);
}




void SampleMaster::printInitializedFiles() const 
{
  bool isRestart = cl.getCheckpointId().compare("") != 0;
  auto initString = isRestart ? "regenerated" : "initialized" ; 

  tout << initString  << " diagnostics file " << diagFile.getFileName()  << std::endl; 

  for(auto &run : runs)
    {
      for(auto &elem : run.getAllFileNames())
	{
	  tout << initString << " file " << elem << std::endl; 
	}
    }
} 


void SampleMaster::printInitialState(const ParallelSetup &pl)
{    
  // compute the initial state 
  for(auto &run: runs)
    {    
      for(nat i = 0; i < run.getChains().size(); ++i)	
	{
	  auto &chain = run.getChains()[i]; 
	  if(pl.isMyChain(run.getRunid(), i))
	    {
	      chain.resume(true, 
#ifdef _DISABLE_INIT_LNL_CHECK
			   false
#else 
			   true 
#endif
			   ); 
	      chain.suspend();
	    }
	}
    }


  // if we are parallel, we must inform the master about what we
  // computed
  pl.synchronizeChainsAtMaster(runs, CommFlag::PrintStat);

  tout << std::endl << "initial state: " << endl; 
  tout << "================================================================" << endl; 
  for(auto &run: runs)
    {
      for(auto &chain: run.getChains())
	{
	  tout << chain ; 
	  tout << "\tRNG(" << chain.getChainRand() << ")" << std::endl; 
	}
      tout << "================================================================" << endl; 
    }
  tout << endl; 

}  


std::pair<double,double> SampleMaster::convergenceDiagnostic(nat &start, nat &end)
{
  if(runParams.getNumRunConv() == 1)
    return std::pair<double,double>(std::numeric_limits<double>::min(),std::numeric_limits<double>::min());

  auto fns = vector<string>{}; 
  for(nat i = 0; i < runParams.getNumRunConv(); ++i)
    {
      auto &&ss = stringstream{}; 
      ss << OutputFile::getFileBaseName(cl.getWorkdir())  << "_topologies." << cl.getRunid() << "." << i; 
      
      // if we do not find the file, let's assume we have multiple branch lengths 
      if (not std::ifstream(ss.str()) ) 
	{
	  ss.str(""); 
	  ss << OutputFile::getFileBaseName(cl.getWorkdir())  << "_topologies." << cl.getRunid() << "." << i << ".tree.0"; 
	}

      fns.push_back(ss.str());
    }
     
  auto&& asdsf = AvgSplitFreqAssessor(fns);

  end = asdsf.getEnd();
      
  int treesInBatch = runParams.getDiagFreq() / runParams.getSamplingFreq(); 

  end /= treesInBatch; 
  end *= treesInBatch;       

  if(end == 0)
    return make_pair(nan(""), nan(""));   

  asdsf.setEnd(end);
  if( runParams.getBurninGen() > 0 )
    {
      assert(runParams.getBurninProportion() == 0.); 

      int treesToDiscard =  runParams.getBurninGen() / runParams.getSamplingFreq(); 

      if(int(end) < treesToDiscard + 2 )
	return make_pair(nan(""), nan("")); 
      else 
	asdsf.setStart(treesToDiscard);  
    }
  else 
    {
      assert(runParams.getBurninGen() == 0); 
      start = (int)((double)end * runParams.getBurninProportion()  ); 
      asdsf.setStart(start);
    } 

  asdsf.extractBipsNew();
  auto asdsfVals = asdsf.computeAsdsfNew(runParams.getAsdsfIgnoreFreq());

  return asdsfVals;
}


void SampleMaster::processConfigFile(string configFileName, const TreeAln* tralnPtr,
				      std::vector<std::unique_ptr<AbstractProposal> > &proposalResult, 
				      std::vector<std::unique_ptr<AbstractParameter> > &variableResult, 
				      std::vector<ProposalSet> &proposalSets,
				      const std::unique_ptr<LikelihoodEvaluator> &eval)
{
  ConfigReader reader; 
  ifstream fh(configFileName); 
  NxsToken token(fh); 

  paramBlock.setTree(tralnPtr); 
  reader.Add(&paramBlock); 

  BlockPrior priorBlock(tralnPtr->getNumberOfPartitions());
  reader.Add(&priorBlock);

  reader.Add(&runParams);

  BlockProposalConfig proposalConfig; 
  reader.Add(&proposalConfig);   

  reader.Execute(token);

  auto r = RunFactory{}; 
  proposalResult = r.produceProposals(proposalConfig, priorBlock , paramBlock, *tralnPtr, eval, runParams.isComponentWiseMH(),  proposalSets);
  variableResult = r.getRandomVariables(); 
}


CLOCK::system_clock::time_point 
SampleMaster::printDuringRun(nat gen,  ParallelSetup &pl)   
{
  pl.synchronizeChainsAtMaster(runs, CommFlag::PrintStat);
  
  stringstream ss ; 
  ss << SOME_FIXED_PRECISION; 
  ss << "[" << gen << "," << CLOCK::duration_cast<CLOCK::duration<double> > (CLOCK::system_clock::now()- lastPrintTime   ).count()     << "s]\t"; 

  bool isFirst = true; 
  for(auto &run : runs)
    {
      if(isFirst)
	isFirst = false; 
      else 
	ss << " ==="; 

      for(auto &c : run.getChains())
	{
	  ss << " " << c.getLikelihood() ; 
	}
    }

  // print it 
  tout << ss.str() << std::endl; 

  return CLOCK::system_clock::now(); 
}

void SampleMaster::run()
{
  tout << "Starting MCMC sampling. Will print the log-likelihoods of\n"
       << "all chains for the given print interval. Independent runs\n" 
       << "are separated by '=', coupled chains are sorted according\n"
       << "to heat (cold chain first)." << std::endl; 
  tout << "First column indicates generation number (completed by all\n"
       << "chains) and the time elapsed while computing these generations "
       << std::endl << std::endl; 
  tout << pl << std::endl; 

  bool hasConverged = false;   
  nat curGen = runs[0].getChains()[0].getGeneration(); // dangerous 
  for(auto & run: runs)
    for(auto &c : run.getChains())
      assert( nat(c.getGeneration()) == curGen );       

  int lastPrint = (curGen / runParams.getPrintFreq() )  * runParams.getPrintFreq() ; 
  int lastDiag =  (curGen / runParams.getDiagFreq() ) * runParams.getDiagFreq(); 
  int lastChkpnt = ( curGen / runParams.getChkpntFreq() )* runParams.getChkpntFreq() ; 
  
  lastPrintTime = printDuringRun(curGen,pl); 

  while(curGen < nat(runParams.getNumGen()) || not hasConverged)   
    { 
      nat nextPrint =  lastPrint + runParams.getPrintFreq(); 
      nat nextDiag = lastDiag + runParams.getDiagFreq(); 
      nat nextChkpnt = lastChkpnt + runParams.getChkpntFreq(); 

      std::vector<nat> stopPoints = { nextChkpnt , nextPrint, nextDiag } ; 

      if(curGen < runParams.getNumGen())
	stopPoints.push_back(runParams.getNumGen());

      nat nextStop = *(std::min_element(stopPoints.begin(), stopPoints.end())); 
      int toExecute = nextStop - curGen; 

      // main part execute 
      for(auto &run : runs)
	{
	  if(pl.isMyRun(run.getRunid()))
	    run.executePart(curGen, toExecute, pl );
	}
      curGen += toExecute; 

      hasConverged = (  runs.size() == 1) && (curGen >= runParams.getNumGen()); 

      if(curGen % runParams.getDiagFreq() == 0 )
	{
	  auto asdsf = make_pair(nan(""), nan("")); 

	  if(runs.size() > 1)
	    {
	      nat start = 0, 
		end = 0; 
	      asdsf = convergenceDiagnostic(start, end); 

	      double convCrit = runParams.getAsdsfConvergence();  
	      if(runParams.isUseAsdsfMax())
		hasConverged = asdsf.second < convCrit;  
	      else 
		hasConverged = asdsf.first < convCrit;  

	      tout << std::endl  << "ASDSF for trees " << start << "-" << end  << " (avg/max):\t"
		   << PERC_PRECISION << asdsf.first * 100 << "%\t" << asdsf.second * 100 << "%"   << std::endl << std::endl; 
	    }

	  pl.synchronizeChainsAtMaster(runs, CommFlag::PrintStat | CommFlag::Swap | CommFlag::Proposals); 
	  if(pl.isGlobalMaster()) 
	    diagFile.printDiagnostics(curGen, asdsf.first, runs);
	  lastDiag = curGen; 
	}
      
      if(curGen % runParams.getPrintFreq() == 0 )
	{
	  lastPrintTime = printDuringRun(curGen,pl); 
	  lastPrint = curGen; 
	}

      if(curGen % runParams.getChkpntFreq() == 0)
	{
	  writeCheckpointMaster();
	  lastChkpnt = curGen; 
	} 
    }

#if defined( _EXPERIMENTAL_INTEGRATION_MODE ) && defined(_GO_TO_INTEGRATION_MODE)

  // go into integration mode, if we want to do so 
  branchLengthsIntegration();
#endif

#if defined(_EXPERIMENTAL_INTEGRATION_MODE)  && defined(_GO_TO_TREE_MOVE_INTEGARTION)
  tInt->integrateAllBranches(*(runs[0].getChains()[0].getTralnPtr()), cl.getRunid());
#endif
  
}
 
 
void SampleMaster::finalizeRuns()
{
  for(auto &run : runs)
    {
      for(auto &chain : run.getChains())
	{
	  if(chain.getChainHeat() == 1. )
	    {
	      tout << "best state was: " << chain.getBestState( )<< endl;       
	    }
	}
      run.finalizeOutputFiles();
    }
  
  tout << endl << "Converged/stopped after " << runs[0].getChains()[0].getGeneration() << " generations" << endl;   
  tout << endl << "Total execution time: " 
       << CLOCK::duration_cast<CLOCK::duration<double> >( CLOCK::system_clock::now() - initTime   ).count() <<  " seconds" << endl; 
}


void SampleMaster::readFromCheckpoint( std::istream &in ) 
{
  for(auto &run : runs)
    run.readFromCheckpoint(in); 

  long start = 0; 
  in >> start; 
  CLOCK::duration<long> durationSinceStart{start};   
  initTime = CLOCK::system_clock::now() -  durationSinceStart; 

}
 
void SampleMaster::writeToCheckpoint( std::ostream &out) const
{  
  for(auto & run : runs)    
    run.writeToCheckpoint(out);

  long duration = CLOCK::duration_cast<CLOCK::duration<long> >(CLOCK::system_clock::now() - initTime).count();
  cWrite(out, duration); 
}



void SampleMaster::writeCheckpointMaster()
{
  pl.synchronizeChainsAtMaster(runs, CommFlag::PrintStat | CommFlag::Proposals | CommFlag::Tree | CommFlag::Swap); 

  if(pl.isGlobalMaster() )
    {

      stringstream ss; 
      ss <<  OutputFile::getFileBaseName(cl.getWorkdir()) << "_newCheckpoint." << cl.getRunid()  ; 
      std::string newName = ss.str();
      ofstream chkpnt; 
      Checkpointable::getOfstream(newName, chkpnt); 
      writeToCheckpoint(chkpnt);
      chkpnt.close(); 

      ss.str("");
      ss << OutputFile::getFileBaseName(cl.getWorkdir()) << "_checkpoint." << cl.getRunid(); 
      std::string curName = ss.str();
      if( std::ifstream(curName) )
	{
	  ss.str("");
	  ss << OutputFile::getFileBaseName(cl.getWorkdir())  << "_prevCheckpointBackup." << cl.getRunid(); 
	  std::string prevName =  ss.str();
	  if( std::ifstream(prevName) )
	    {
	      int ret = remove(prevName.c_str()); 
	      assert(ret == 0); 
	    }
	      
	  int ret = rename(curName.c_str(), prevName.c_str()); 
	  assert(ret == 0); 
	}
	  
      int ret = rename(newName.c_str(), curName.c_str()); 
      assert(ret == 0); 
    }
  else 
    {
      std::ofstream nullstream("/dev/null"); 
      writeToCheckpoint(nullstream); 
      nullstream.close(); 
    }
} 


 #include "IntegrationModuleImpl.hpp"
