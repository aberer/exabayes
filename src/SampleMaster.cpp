#include <sstream>
#include <fstream>
#include <memory>

#include "TreeResource.hpp"
#include "ByteFileResource.hpp"

#include "AdHocIntegrator.hpp"
#include "ProposalSet.hpp"

#include "tree-parse/BasicTreeReader.hpp"
#include "common.h"

#include "config/BlockProposalConfig.hpp"
#include "config/BlockRunParameters.hpp"
#include "config/ConfigReader.hpp"

#include "parameters/TopologyParameter.hpp"

#include "TreeInitializer.hpp"

#include "config/MemoryMode.hpp"
#include "SampleMaster.hpp"
#include "Chain.hpp"
#include "TreeRandomizer.hpp"
#include "tune.h"
#include "RunFactory.hpp"	

#include "GlobalVariables.hpp"

#include "BoundsChecker.hpp"
#include "eval/ArrayRestorer.hpp"
#include "SplitFreqAssessor.hpp"

#include "eval/FullCachePolicy.hpp"
#include "eval/NoCachePolicy.hpp"

#include "proposers/AbstractProposer.hpp"
#include "proposers/SlidingProposal.hpp"
#include "proposers/MultiplierProposal.hpp"

#include "TreeIntegrator.hpp"

#include "parser/PhylipParser.hpp"

// a developmental mode to integrate over branch lengths
#define _GO_TO_INTEGRATION_MODE


#ifndef _EXPERIMENTAL_INTEGRATION_MODE
#undef _GO_TO_INTEGRATION_MODE
#endif

void genericExit(int code); 


SampleMaster::SampleMaster(const ParallelSetup &pl, const CommandLine& cl ) 
  : _pl(pl)
  , _initTime(CLOCK::system_clock::now())
  , _cl(cl)
  , _lastPrintTime(CLOCK::system_clock::now())
{
}


bool SampleMaster::initializeTree(TreeAln &traln, std::string startingTree, Randomness &treeRandomness, const std::vector<AbstractParameter*> &params)
{  
  bool hasBranchLength = false; 
  if(startingTree.compare("") != 0 )
    {
      hasBranchLength = std::any_of(startingTree.begin(), startingTree.end(), [](const char c ){ return c == ':'; }  ); 

      auto &&iss = std::istringstream {startingTree };
      auto reader = BasicTreeReader<NameLabelReader,ReadBranchLength>{traln.getNumberOfTaxa()};

      auto mapAsVect = traln.getTaxa(); 
      mapAsVect.erase(mapAsVect.begin()); 
      auto map = std::unordered_map<std::string,nat>{}; 
      nat ctr = 1; 
      for(auto elem : mapAsVect)
	{
	  map[elem] = ctr; 
	  ++ctr; 
	}
      reader.setLabelMap(map);
      auto branches = reader.extractBranches(iss);
      assert(branches.size() == traln.getNumberOfBranches() ); 

      traln.unlinkTree();
      for(auto b : branches)
	{
	  traln.clipNode(traln.getUnhookedNode(b.getPrimNode()) , traln.getUnhookedNode(b.getSecNode())); 
	  if(hasBranchLength)
	    {
	      for(auto param : params)
		{
		  auto bCopy = b; 
		  bCopy.setConvertedInternalLength(traln, param, b.getLength()); 
		  traln.setBranch(bCopy, param); 
		}
	    }
	}
    }
  else
    {	      
      if(_runParams.isUseParsimonyStarting())
	TreeRandomizer::createParsimonyTree(traln, treeRandomness); 
      else
	TreeRandomizer::randomizeTree(traln, treeRandomness); 
    }

  return hasBranchLength; 
}



std::vector<bool> SampleMaster::initTrees(vector<shared_ptr<TreeAln> > &trees, randCtr_t seed, std::vector<std::string> startingTreeStrings, const std::vector<AbstractParameter*> &params)
{  
  auto hasBl = std::vector<bool>{}; 

  nat treesConsumed = 0; 
  auto treeRandomness = Randomness(seed); 

  auto treesToInitialize = vector<shared_ptr<TreeAln> >{}; 
  treesToInitialize.push_back(trees[0]); 
  if(not _runParams.isHeatedChainsUseSame())
    {
      for(auto iter =  trees.begin() + 1  ; iter < trees.end(); ++iter)
      treesToInitialize.push_back(*iter); 
    }

  // choose how to initialize the topology
  for(auto &tralnPtr : treesToInitialize)
    {
      auto stTr = treesConsumed < startingTreeStrings.size() ? startingTreeStrings.at(treesConsumed) : std::string{""}; 
      ++treesConsumed;
      
      auto hadBl = initializeTree(*tralnPtr, stTr, treeRandomness, params); 
      hasBl.push_back(hadBl); 
    }

  // propagate the tree to the coupled chains, if necessary
  if(_runParams.isHeatedChainsUseSame())
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
  
  return hasBl; 
}


void SampleMaster::printAlignmentInfo(const TreeAln &traln)
{
  tout << std::endl << "The (binary) alignment file you provided, consists of the following\n" 
       << "partitions:" 
       << std::endl;
  
  for(nat i = 0 ;i < traln.getNumberOfPartitions() ;++i)
    {
      auto& partition = traln.getPartition(i);       
      nat length = partition.upper - partition.lower; 
      
      tout << std::endl; 

      tout << "number:\t\t" << i << std::endl;     
      tout << "name:\t\t" << partition.partitionName << std::endl; 
      tout << "#patterns:\t" << length << std::endl;       

      switch(partition.dataType)
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
  if(_pl.getRunsParallel()  > 1 )
    {
      tout << "Will execute "<< _pl.getRunsParallel() << " runs in parallel." << std::endl;       
      if(nat (_runParams.getNumRunConv()) < _pl.getRunsParallel())
	{
	  tout 
	    << "Error: in the configuration file you specified to run " <<  _runParams.getNumRunConv() << " independent\n" 
	    << "runs. Your command line indicates, that you want to run " << _pl.getRunsParallel() << " of them\n"
	    << "in parallel. To shield you from surprises, " << PROGRAM_NAME << " will conservatively abort."<< std::endl; 
	  ParallelSetup::genericExit(-1); 
	}
    }
  if(_pl.getChainsParallel() > 1)
    {
      tout << "Will execute " << _pl.getChainsParallel() << " chains in parallel."<< std::endl; 
      if(nat(_runParams.getNumCoupledChains()) < _pl.getChainsParallel())	
	{
	  tout
	    << "Error: in the configuration file you specified to run " <<  _runParams.getNumCoupledChains() << " coupled\n" 
	    << "chains per run. Your command line indicates, that you want to run " << _pl.getChainsParallel() << " of them\n"
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
  if(_cl.getCheckpointId().compare("") != 0)
    {
      auto prevId = _cl.getCheckpointId(); 
      
      std::stringstream ss ;       
      ss << _cl.getWorkdir() << ( _cl.getWorkdir().compare("") == 0  ? "" : "/") 
	 << PROGRAM_NAME << "_checkpoint." << _cl.getCheckpointId();       
      auto checkPointFile = ss.str(); 
      std::ifstream chkpnt; 
      Serializable::getIfstream(checkPointFile, chkpnt); 
      if( not chkpnt )
	{
	  tout  << "Warning! Could not open / find file " << checkPointFile << std::endl; 
	  tout << "Although extremely unlikely, the previous run may have been killed\n"
	       << "during the writing process of the checkpoint. Will try to recover from backup..." << std::endl; 
	  
	  ss.str(""); 
	  ss << PROGRAM_NAME << "_prevCheckpointBackup." << _cl.getCheckpointId();
	  checkPointFile = ss.str();	  
	  Serializable::getIfstream(checkPointFile, chkpnt); 
	  
	  if( not chkpnt)
	    {
	      tout << "Could not recover checkpoint file backup either. Probably there is no backup, since not enough generations have completed. Giving up. Please start a new run." << std::endl; 
	      ParallelSetup::genericExit(-1); 
	    }
	  else 
	    tout << "Success! You lost one checkpoint, but we can continue from the backup. " << std::endl; 	  
	}

      deserialize(chkpnt); 
      if(_pl.isGlobalMaster())
	{
	  for(auto &run : _runs)
	    run.regenerateOutputFiles(_cl.getWorkdir(), _cl.getCheckpointId());
	  nat curGen = _runs[0].getChains()[0].getGeneration();
	  _diagFile.regenerate(  _cl.getWorkdir(), _cl.getRunid(),  _cl.getCheckpointId(), 
				curGen, _runs[0].getSwapInfo().getMatrix().size() * _runs.size());
	}
    }
}


// TODO could move this method somewhere else  
LikelihoodEvaluator
SampleMaster::createEvaluatorPrototype(const TreeAln &initTree, std::string binaryFile )
{
  auto &&plcy =  std::unique_ptr<ArrayPolicy>();

  switch(_cl.getMemoryMode())
    {
    case MemoryMode::RESTORE_ALL: 
      {
	plcy = std::unique_ptr<ArrayPolicy>(new FullCachePolicy(initTree, true, true));
	if(  (initTree.getMode() & RunModes::MEMORY_SEV) != RunModes::NOTHING )
	  plcy->enableRestoreGapVector();
      }
      break; 
    case MemoryMode::RESTORE_INNER_TIP: 
      {
	plcy = std::unique_ptr<ArrayPolicy>(new FullCachePolicy(initTree, false, true));
	if(  (initTree.getMode() & RunModes::MEMORY_SEV) != RunModes::NOTHING )
	  plcy->enableRestoreGapVector();
      }
      break; 
    case MemoryMode::RESTORE_INNER_INNER: 
      {
	plcy = std::unique_ptr<ArrayPolicy>(new FullCachePolicy(initTree, false, false));
	if(  (initTree.getMode() & RunModes::MEMORY_SEV) != RunModes::NOTHING )
	  plcy->enableRestoreGapVector();
      }
      break; 
    case MemoryMode::RESTORE_NONE: 
      {
	plcy = std::unique_ptr<ArrayPolicy>(new NoCachePolicy(initTree)); 
      }
      break; 
    default: 
      assert(0); 
    }

  auto eval = LikelihoodEvaluator(initTree, plcy.get()); 

#ifdef DEBUG_LNL_VERIFY

  auto dT = make_shared<TreeAln>(initTree.getNumberOfTaxa());
  auto &&ti = TreeInitializer(std::unique_ptr<InitializationResource>(new ByteFileResource(binaryFile))); 
  ti.initializeWithAlignmentInfo(*dT, RunModes::NOTHING); 



  // dT->initializeFromByteFile(binaryFile, RunModes::NOTHING ); 

  // dT->enableParsimony(); 
  eval.setDebugTraln(dT);
#endif

  return eval; 
}


void SampleMaster::initializeWithParamInitValues(std::vector<shared_ptr<TreeAln>> &tralns , 
						 const std::vector<AbstractParameter*> &params,
						 const std::vector<bool> hasBls) const 
{
  //  do normal parameters first 
  for(auto &param : params)
    {
      auto cat = param->getCategory( ); 
      if( cat != Category::TOPOLOGY && cat != Category::BRANCH_LENGTHS)
	{

	  // :NOTICE: treat prot frequencies differently!
	  if(cat == Category::FREQUENCIES
	     && tralns[0]->getPartition(param->getPartitions()[0] ) .dataType == AA_DATA ) 
	    {
	      for(auto p : param->getPartitions())
		{
		  for(auto &tralnPtr : tralns ) 
		    {
		      auto& partition = tralnPtr->getPartition(p);
		      partition.protFreqs = TRUE; 
		    }
		}
	    }
	  
	  // :NOTICE: rev mat parameters for aa-partitions must be prepared here! 
	  if(cat == Category::SUBSTITUTION_RATES
	     && tralns[0]->getPartition(param->getPartitions()[0]).dataType == AA_DATA)
	    {
	      for(auto p : param->getPartitions())
		{
		  for(auto &tralnPtr : tralns)
		    {
		      auto &partition = tralnPtr->getPartition(p);
		      partition.protModels = GTR; 	
		    }
		}
	    }

	  auto&& prior = param->getPrior(); 
	  auto content = prior->getInitialValue();
	  
	  auto &bla = *(tralns[0].get()); 
	  param->verifyContent(bla, content); 

	  for(auto &traln : tralns)
	    param->applyParameter(*traln, content); 

	  auto result = param->extractParameter(*(tralns[0])); 
	}
    }

  for(auto &param : params)
    {
      auto cat = param->getCategory(); 
      if(cat == Category::BRANCH_LENGTHS)
	{
	  auto &&prior = param->getPrior();
	  auto content = prior->getInitialValue();

	  auto ctr = nat{0};
	  for(auto &tralnPtr : tralns)
	    {
	      if( not hasBls[ctr])
		{
		  for(auto &b : tralnPtr->extractBranches(param))
		    {
		      b.setConvertedInternalLength( *tralnPtr,param,  content.values[0] );

		      if(not BoundsChecker::checkBranch(b))
			BoundsChecker::correctBranch(b); 

		      tralnPtr->setBranch(b, param);
		    }
		}
	    }
	}
    }
}



std::vector<std::string> SampleMaster::getStartingTreeStrings()
{
  auto result =  std::vector<std::string>{};

  auto&& ifh = std::ifstream{_cl.getTreeFile()}; 
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
  
  tout << "Will employ the following proposal mixture (frequency,id,type,affected variables): " << endl; 
  for(auto &p : proposals )
    {
      tout << PERC_PRECISION << p->getRelativeWeight() / sum * 100 <<   "%\t" ; 
      tout << p->getId() << "\t" ; 
      p->printShort(tout ) ; 
      tout << endl; 
    }
  if(proposals.size() == 0)
    tout << "None." << std::endl; 

  tout << std::endl; 

  if(_runParams.isComponentWiseMH())
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

      tout << "\tinit value:\t" ; 

      if(dynamic_cast<TopologyParameter*>(v.get())  != nullptr )
	{
	  auto p = v->getPrior(); 

	  if(not p->needsIntegration())
	    tout << "fixed" ; 
	  else if(_runParams.isUseParsimonyStarting())
	    tout << "parsimony" ; 
	  else 
	    tout << "random" ; 
	  
	  if(_cl.getTreeFile().compare("") != 0)
	    tout << " or given (tree file)" ; 
	  tout << std::endl; 
	}
      else 
	tout << v->getPrior()->getInitialValue() << std::endl; 
    }
  tout << "================================================================" << std::endl;
  tout << endl; 
}


std::string SampleMaster::getOrCreateBinaryFile() const 
{
  auto binaryAlnFile =std::string{}; 
  if( not _cl.alnFileIsBinary())
    {
      tout << "You provided an alignment file in phylip format. Trying to parse it..." << std::endl; 
      
      bool haveModelFile = _cl.getModelFile().compare("") != 0; 
      auto modelInfo = haveModelFile ? _cl.getModelFile( ): _cl.getSingleModel(); 

      auto parser = PhylipParser{ _cl.getAlnFileName() , modelInfo, haveModelFile}; 
      
      binaryAlnFile = std::string(_cl.getWorkdir() 
					+  ( _cl.getWorkdir().compare("") == 0 ? "" : "/"   ) 
					+   "ExaBayes_binaryAlignment" + "." + _cl.getRunid()) ;
      if(std::ifstream(binaryAlnFile))
	{
	  tout << "removing previous binary alignment representation " << binaryAlnFile << std::endl; 
	  remove(std::string(binaryAlnFile).c_str()); 
	}
      
      parser.parse(); 
      parser.writeToFile(binaryAlnFile); 
      // tout << "wrote to "<< binaryAlnFile << std::endl; 
      
      if(not std::ifstream(binaryAlnFile))
	{
	  tout << "Error: tried to create intermediate file "<< binaryAlnFile << ", but did not succeed!"  << std::endl; 
	  exit(-1); 
	}
    }
  else 
    {
      binaryAlnFile = _cl.getAlnFileName();   
    }
  
  return binaryAlnFile; 
}



nat SampleMaster::peekNumTax(std::string filePath)
{
  nat result = 0; 
  auto &&in = std::ifstream{filePath, std::ios::binary}; 
  auto initId = std::string{"BINARY"}; 
  
  for(nat i = 0; i < initId.size(); ++i)
    {
      int ch = 0; 
      ch = in.get(); 

      if(ch != int(initId[i]))
	{
	  tout << "error: expected >" << char(initId[i])  << "< but got >" << ch << "<"  << std::endl; 
	  exit(-1); 
	}
    }

  in.read((char*)&result, sizeof(result)) ; 

  return result ;
}


void SampleMaster::initializeRuns(Randomness rand)
{  
  auto startingTrees = getStartingTreeStrings(); 

  if(_cl.getCheckpointId().compare(_cl.getRunid()) == 0)
    {
      std::cerr << "You specified >" << _cl.getRunid() << "< as runid and intended\n"
		<< "to restart from a previous run with id >" << _cl.getCheckpointId() << "<."
		<< "Please specify a new runid for the restart. " << std::endl; 
      ParallelSetup::genericExit(-1); 
    }

  auto binaryAlnFile = getOrCreateBinaryFile(); 
  auto numTax = peekNumTax(binaryAlnFile); 

  // initialize one tree 
  auto&& initTreePtr = std::shared_ptr<TreeAln>(new TreeAln(numTax)); 

  auto runmodes = _cl.getTreeInitRunMode();

  auto trees =  std::vector<std::shared_ptr<TreeAln> >{}; 
  
  auto &&ti = TreeInitializer(std::unique_ptr<InitializationResource>(new ByteFileResource(binaryAlnFile))); 
  ti.initializeWithAlignmentInfo(*initTreePtr, runmodes); 

  trees.push_back(initTreePtr);
  auto &initTree = *(trees[0].get()); 

  // START integrator
#ifdef _EXPERIMENTAL_INTEGRATION_MODE
  std::shared_ptr<TreeAln> aTree = std::unique_ptr<TreeAln>(new TreeAln()); 
  aTree->initializeFromByteFile(binaryAlnFile, runmodes); 
  aTree->enableParsimony();

  // let's have another tree for debug
  auto dT = make_shared<TreeAln>();
  dT->initializeFromByteFile(binaryAlnFile, runmodes); 
  dT->enableParsimony(); 

  TreeRandomizer::randomizeTree(*aTree, masterRand); 
  ahInt = new AdHocIntegrator(aTree, dT, masterRand.generateSeed());

  tInt = new TreeIntegrator(aTree, dT, masterRand.generateSeed()); 
#endif
  // END

 auto evalUptr = createEvaluatorPrototype(initTree,  binaryAlnFile); 
  
  auto&& proposals =  std::vector<unique_ptr<AbstractProposal> >{} ; 
  auto&& params = std::vector<unique_ptr<AbstractParameter> >{} ; 
  auto&& proposalSets =  std::vector<ProposalSet>{};  
  std::tie(params, proposals, proposalSets) = processConfigFile(_cl.getConfigFileName(), initTree);

  assert(_runParams.getTuneFreq() > 0); 

  // ORDER: must be after initWithConfigFile

  for(nat i = 1 ; i < _runParams.getNumCoupledChains(); ++i)
    {      
      trees.push_back(std::make_shared<TreeAln>(numTax)); 

      tout << "\n\n copying a tree\n\n " ; 

      auto &&ti = TreeInitializer(std::unique_ptr<InitializationResource>(new TreeResource(initTreePtr.get())));
      ti.initializeWithAlignmentInfo(*(trees[i]), runmodes); 

#if HAVE_PLL == 0
      if(_cl.isPerPartitionDataDistribution())
	{
	  auto &tr = trees[i]->getTrHandle(); 
	  tr.manyPartitions = TRUE; 
	}
#endif
    }

  auto runSeeds = vector<randCtr_t>{};
  auto treeSeeds = vector<randCtr_t>{}; 
  for(nat i = 0; i < _runParams.getNumRunConv();++i)
    {
      runSeeds.push_back(rand.generateSeed()); 
      treeSeeds.push_back(rand.generateSeed()); 
    }

  // determine if topology is fixed 
  bool topoIsFixed = false; 
  for(auto &v :params)
    {
      if( dynamic_cast<TopologyParameter*>(v.get()) != nullptr
	  && not v->getPrior()->needsIntegration() )
	topoIsFixed = true; 
    }
  
  // gather branch length parameters
  auto blParams = std::vector<AbstractParameter*>{} ;
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
      auto treeRandomness = Randomness(treeSeeds[0]); 
      auto &something = *initTreePtr; 
      initializeTree(something, 
		     startingTrees.size() > 0  ? startingTrees.at(0): std::string(""), // meh 
		     treeRandomness, blParams); 
    }

  for(nat i = 0; i < _runParams.getNumRunConv() ; ++i)
    {    
      auto hadBls = std::vector<bool>(trees.size() , false);
      if(topoIsFixed)
	{
	  for(auto &t : trees)
	    t->copyModel(*initTreePtr); 
	}
      else 	
	hadBls = initTrees(trees, treeSeeds[i],  startingTrees, blParams); 

      auto paramView = std::vector<AbstractParameter*>{}; 
      for(auto &param : params)
	paramView.push_back(param.get()); 

      initializeWithParamInitValues(trees, paramView, hadBls );
      
      auto chains = vector<Chain>{};       
      for(nat j = 0; j < _runParams.getNumCoupledChains(); ++j)
	{
	  randCtr_t c; 
	  auto &t = trees.at(j);
	  chains.emplace_back( c , t, proposals, proposalSets, evalUptr, _cl.isDryRun() ); 
	  auto &chain = chains[j]; 		
	  chain.setRunId(i); 
	  chain.setTuneFreuqency(_runParams.getTuneFreq()); 
	  chain.setHeatIncrement(j); 
	  chain.setDeltaT(_runParams.getHeatFactor()); 
	}
	  
      _runs.emplace_back(runSeeds[i], i, _cl.getWorkdir(), _cl.getRunid(), _runParams.getNumCoupledChains(), chains); 
      auto &run = *(_runs.rbegin()); 
      run.setTemperature(_runParams.getHeatFactor());
      run.setPrintFreq(_runParams.getPrintFreq()); 
      run.setSwapInterval(_runParams.getSwapInterval()); 
      run.setSamplingFreq(_runParams.getSamplingFreq()); 
      run.setNumSwaps(_runParams.getNumSwaps());
      run.seedChains(); 	  
      
      if(_pl.isRunLeader() && _pl.isMyRun(run.getRunid()))
	run.initializeOutputFiles(_cl.isDryRun());
    }

  initializeFromCheckpoint(); 

  if(_pl.isGlobalMaster() && not _diagFile.isInitialized() && not _cl.isDryRun() )
    _diagFile.initialize(_cl.getWorkdir(), _cl.getRunid(), _runs);

  // post-pone all printing to the end  
  printAlignmentInfo(initTree);   
  printParameters(initTree, params);
  printProposals(proposals, proposalSets); 
  informPrint();
  
  printInitializedFiles(); 

  if(not _cl.isDryRun())
    printInitialState(); 
  
  _pl.printLoadBalance(initTree);
}


void SampleMaster::printInitializedFiles() const 
{
  bool isRestart = _cl.getCheckpointId().compare("") != 0;
  auto initString = isRestart ? "regenerated" : "initialized" ; 

  tout << initString  << " diagnostics file " << _diagFile.getFileName()  << std::endl; 

  for(auto &run : _runs)
    {
      for(auto &elem : run.getAllFileNames())
	{
	  tout << initString << " file " << elem << std::endl; 
	}
    }
} 


void SampleMaster::printInitialState()  
{    
  // compute the initial state 
  for(auto &run: _runs)
    {    
      for(nat i = 0; i < run.getChains().size(); ++i)	
	{
	  auto &chain = run.getChains()[i]; 
	  if(_pl.isMyChain(run.getRunid(), i))
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
  _pl.synchronizeChainsAtMaster(_runs, CommFlag::PrintStat);

  tout << std::endl << "initial state: " << endl; 
  tout << "================================================================" << endl; 
  for(auto &run: _runs)
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
  if(_runParams.getNumRunConv() == 1)
    return std::pair<double,double>(std::numeric_limits<double>::min(),std::numeric_limits<double>::min());

  auto fns = vector<string>{}; 
  for(nat i = 0; i < _runParams.getNumRunConv(); ++i)
    {
      auto &&ss = stringstream{}; 
      ss << OutputFile::getFileBaseName(_cl.getWorkdir())  << "_topologies." << _cl.getRunid() << "." << i; 
      
      // if we do not find the file, let's assume we have multiple branch lengths 
      if (not std::ifstream(ss.str()) ) 
	{
	  ss.str(""); 
	  ss << OutputFile::getFileBaseName(_cl.getWorkdir())  << "_topologies." << _cl.getRunid() << "." << i << ".tree.0"; 
	}

      fns.push_back(ss.str());
    }
     
  auto&& asdsf = SplitFreqAssessor(fns);

  end = asdsf.getMinNumTrees();   

  int treesInBatch = _runParams.getDiagFreq() / _runParams.getSamplingFreq(); 
  end /= treesInBatch; 
  end *= treesInBatch;       

  if(end == 0)
    return make_pair(nan(""), nan(""));   

  if( _runParams.getBurninGen() > 0 )
    {
      assert(_runParams.getBurninProportion() == 0.); 

      int treesToDiscard =  _runParams.getBurninGen() / _runParams.getSamplingFreq(); 

      if(int(end) < treesToDiscard + 2 )
	return make_pair(nan(""), nan("")); 
      else 
	start = treesToDiscard; 
    }
  else 
    {
      assert(_runParams.getBurninGen() == 0); 
      start = (int)((double)end * _runParams.getBurninProportion()  ); 
    } 

  asdsf.extractBipsNew(start, end, false);
  auto asdsfVals = asdsf.computeAsdsfNew(_runParams.getAsdsfIgnoreFreq());

  return asdsfVals;
}


auto SampleMaster::processConfigFile(string configFileName, const TreeAln &traln )
  ->   std::tuple<std::vector<std::unique_ptr<AbstractParameter> > , std::vector<std::unique_ptr<AbstractProposal> > , std::vector<ProposalSet> >  
{
  auto reader = ConfigReader{}; 
  auto && fh = ifstream(configFileName); 
  auto token = NxsToken(fh); 
  auto&& paramBlock = BlockParams{}; 
  paramBlock.setTree(&traln); 
  auto&& priorBlock = BlockPrior(traln.getNumberOfPartitions());
  auto proposalConfig = BlockProposalConfig{}; 
  
  reader.Add(&paramBlock); 
  reader.Add(&priorBlock);
  reader.Add(&_runParams);
  reader.Add(&proposalConfig);   
  reader.Execute(token);

  auto r = RunFactory{}; 
  auto paramResult = paramBlock.getParameters();

  r.addStandardParameters(paramResult, traln);
  auto proposalSetResult = std::vector<ProposalSet>{};
  auto&& proposalResult = std::vector<std::unique_ptr<AbstractProposal>>{}; 
  std::tie(proposalResult, proposalSetResult) = r.produceProposals(proposalConfig, priorBlock , paramResult, traln, _runParams.isComponentWiseMH());

  // sanity check 
  for(auto &paramPtr : paramResult)
    paramPtr->checkSanityPartitionsAndPrior(traln);

  // now enumerate the proposals 
  nat ctr = 0; 
  for(auto &p : proposalResult)
    {
      p->setId(ctr); 
      ++ctr; 
    }
  for(auto &p : proposalSetResult)
    ctr = p.numerateProposals(ctr);

  return std::make_tuple(std::move(paramResult), std::move(proposalResult), proposalSetResult);
}


CLOCK::system_clock::time_point 
SampleMaster::printDuringRun(nat gen)   
{
  _pl.synchronizeChainsAtMaster(_runs, CommFlag::PrintStat);
  
  stringstream ss ; 
  ss << SOME_FIXED_PRECISION; 
  ss << "[" << gen << "," << CLOCK::duration_cast<CLOCK::duration<double> > (CLOCK::system_clock::now()- _lastPrintTime   ).count()     << "s]\t"; 

  bool isFirst = true; 
  for(auto &run : _runs)
    {
      if(isFirst)
	isFirst = false; 
      else 
	ss << " ==="; 
      
      auto sortedLnls = std::vector<std::pair<nat,double>>{}; 
      for(auto &c : run.getChains() ) 
	sortedLnls.emplace_back(c.getCouplingId(), c.getLikelihood()); 
      std::sort(sortedLnls.begin(), sortedLnls.end(), [] (const std::pair<nat,double> &elem1, const std::pair<nat,double> &elem2 ) { return elem1.first < elem2.first;  }); 

      for(auto elem : sortedLnls)
	ss << " " << elem.second; 
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
  tout << _pl << std::endl; 

  bool hasConverged = false;   
  nat curGen = _runs[0].getChains()[0].getGeneration(); // dangerous 
  for(auto & run: _runs)
    for(auto &c : run.getChains())
      assert( nat(c.getGeneration()) == curGen );       

  int lastPrint = (curGen / _runParams.getPrintFreq() )  * _runParams.getPrintFreq() ; 
  int lastDiag =  (curGen / _runParams.getDiagFreq() ) * _runParams.getDiagFreq(); 
  int lastChkpnt = ( curGen / _runParams.getChkpntFreq() )* _runParams.getChkpntFreq() ; 
  
  _lastPrintTime = printDuringRun(curGen); 

  while(curGen < nat(_runParams.getNumGen()) || not hasConverged)   
    { 
      nat nextPrint =  lastPrint + _runParams.getPrintFreq(); 
      nat nextDiag = lastDiag + _runParams.getDiagFreq(); 
      nat nextChkpnt = lastChkpnt + _runParams.getChkpntFreq(); 

      std::vector<nat> stopPoints = { nextChkpnt , nextPrint, nextDiag } ; 

      if(curGen < _runParams.getNumGen())
	stopPoints.push_back(_runParams.getNumGen());

      nat nextStop = *(std::min_element(stopPoints.begin(), stopPoints.end())); 
      int toExecute = nextStop - curGen; 

      // main part execute 
      for(auto &run : _runs)
	{
	  if(_pl.isMyRun(run.getRunid()))
	    run.executePart(curGen, toExecute, _pl );
	}
      curGen += toExecute; 

      hasConverged = (  _runs.size() == 1) && (curGen >= _runParams.getNumGen()); 

      if(curGen % _runParams.getDiagFreq() == 0 )
	{
	  auto asdsf = make_pair(nan(""), nan("")); 

	  if(_runs.size() > 1)
	    {
	      nat start = 0, 
		end = 0; 
	      asdsf = convergenceDiagnostic(start, end); 

	      double convCrit = _runParams.getAsdsfConvergence();  
	      if(_runParams.isUseAsdsfMax())
		hasConverged = asdsf.second < convCrit;  
	      else 
		hasConverged = asdsf.first < convCrit;  

	      tout << std::endl  << "ASDSF for trees " << start << "-" << end  << " (avg/max):\t"
		   << PERC_PRECISION << asdsf.first * 100 << "%\t" << asdsf.second * 100 << "%"   << std::endl << std::endl; 
	    }

	  _pl.synchronizeChainsAtMaster(_runs, CommFlag::PrintStat | CommFlag::Swap | CommFlag::Proposals); 
	  if(_pl.isGlobalMaster()) 
	    _diagFile.printDiagnostics(curGen, asdsf.first, _runs);
	  lastDiag = curGen; 
	}
      
      if(curGen % _runParams.getPrintFreq() == 0 )
	{
	  _lastPrintTime = printDuringRun(curGen); 
	  lastPrint = curGen; 
	}

      if(curGen % _runParams.getChkpntFreq() == 0)
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
  auto sprLengthStr = std::string(std::getenv("SPR_LENGTH")); 
  auto && iss  = std::istringstream{sprLengthStr}; 
  nat sprLen = 0;  
  iss >> sprLen; 
  
  tInt->integrateAllBranchesNew(*(_runs[0].getChains()[0].getTralnPtr()), 
				_cl.getRunid(),
				sprLen );
#endif
}

 
void SampleMaster::finalizeRuns()
{
  for(auto &run : _runs)
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
  
  tout << endl << "Converged/stopped after " << _runs[0].getChains()[0].getGeneration() << " generations" << endl;   
  tout << endl << "Total execution time: " 
       << CLOCK::duration_cast<CLOCK::duration<double> >( CLOCK::system_clock::now() - _initTime   ).count() <<  " seconds" << endl; 
}


void SampleMaster::deserialize( std::istream &in ) 
{
  for(auto &run : _runs)
    run.deserialize(in); 

  long start = 0; 
  in >> start; 
  CLOCK::duration<long> durationSinceStart{start};   
  _initTime = CLOCK::system_clock::now() -  durationSinceStart; 

}
 
void SampleMaster::serialize( std::ostream &out) const
{  
  for(auto & run : _runs)    
    run.serialize(out);

  long duration = CLOCK::duration_cast<CLOCK::duration<long> >(CLOCK::system_clock::now() - _initTime).count();
  cWrite(out, duration); 
}



void SampleMaster::writeCheckpointMaster()
{
  _pl.synchronizeChainsAtMaster(_runs, CommFlag::PrintStat | CommFlag::Proposals | CommFlag::Tree | CommFlag::Swap); 

  if(_pl.isGlobalMaster() )
    {

      stringstream ss; 
      ss <<  OutputFile::getFileBaseName(_cl.getWorkdir()) << "_newCheckpoint." << _cl.getRunid()  ; 
      std::string newName = ss.str();
      ofstream chkpnt; 
      Serializable::getOfstream(newName, chkpnt); 
      serialize(chkpnt);
      chkpnt.close(); 

      ss.str("");
      ss << OutputFile::getFileBaseName(_cl.getWorkdir()) << "_checkpoint." << _cl.getRunid(); 
      std::string curName = ss.str();
      if( std::ifstream(curName) )
	{
	  ss.str("");
	  ss << OutputFile::getFileBaseName(_cl.getWorkdir())  << "_prevCheckpointBackup." << _cl.getRunid(); 
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
      auto &&nullstream = std::ofstream("/dev/null"); 
      serialize(nullstream); 
      nullstream.close(); 
    }
} 


#include "IntegrationModuleImpl.hpp"
