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


// a developmental mode to integrate over branch lengths
// #define _GO_TO_INTEGRATION_MODE

void genericExit(int code); 
static int countNumberOfTreesQuick(const char *fn ); 

SampleMaster::SampleMaster(const ParallelSetup &pl, const CommandLine& _cl ) 
  : pl(pl)
  , initTime(CLOCK::system_clock::now())
  , masterRand(cl.getSeed())
  , cl(_cl)
  , timeIncrement(CLOCK::system_clock::now())
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
	  tralnPtr->setBranch(b); 
	}
    }

  // propagate the tree to the coupled chains, if necessaroy
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



void SampleMaster::printAlignmentInfo(const TreeAln &traln)
{
  tout << std::endl << "The (binary) alignment file you provided, consists of the following\n" 
       << "partitions:" 
       << std::endl;
  
  for(int i = 0 ;i < traln.getNumberOfPartitions() ;++i)
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

  tout << "\nParameter names will reference the above partition numbers."<< std::endl  ; 
}


void SampleMaster::initializeRuns( )
{  
  FILE *treeFH = NULL; 
  nat numTreesAvailable = countNumberOfTreesQuick(cl.getTreeFile().c_str()); 

  if(cl.getCheckpointId().compare(cl.getRunid()) == 0)
    {
      std::cerr << "You specified >" << cl.getRunid() << "< as runid and inteded to restart from a previous run with id >" << cl.getCheckpointId() << "<. Please specify a new runid for the restart. " << std::endl; 
      ParallelSetup::genericExit(-1); 
    }

  // initialize one tree 
  vector<shared_ptr<TreeAln> > trees; 
  trees.push_back(make_shared<TreeAln>());
  trees[0]->initializeFromByteFile(cl.getAlnFileName()); 
  trees[0]->enableParsimony();

  printAlignmentInfo(*(trees[0])); 

  vector<unique_ptr<AbstractProposal> > proposals; 
  vector<unique_ptr<AbstractParameter> > variables; 

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
	      chains.emplace_back( c , trees[j], proposals, eval  ); 
	      auto &chain = chains[j]; 		
	      chain.setRunId(i); 
	      chain.setTuneFreuqency(runParams.getTuneFreq()); 
	      chain.setHeatIncrement(j); 
	      chain.setDeltaT(runParams.getHeatFactor()); 
	    }
	  
	  runs.emplace_back(runSeeds[i], i, cl.getWorkdir(), runParams.getNumCoupledChains(), chains); 
	  auto &run = *(runs.rbegin()); 
	  run.setTemperature(runParams.getHeatFactor());
	  run.setTuneHeat(runParams.getTuneHeat()); 
	  run.setPrintFreq(runParams.getPrintFreq()); 
	  run.setSwapInterval(runParams.getSwapInterval()); 
	  run.setSamplingFreq(runParams.getSamplingFreq()); 
	  run.setRunName(cl.getRunid()); 
	  run.seedChains(); 	  
	  
	  if(pl.isReportingProcess())
	    run.initializeOutputFiles ();
	}
    }
  
  // continue from checkpoint  
  if(cl.getCheckpointId().compare("") != 0)
    {
      auto prevId = cl.getCheckpointId(); 

      std::stringstream ss ; 
      ss << PROGRAM_NAME << "_checkpoint." << cl.getCheckpointId();       
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
      if(pl.isMasterReporter())
	{
	  for(auto &run : runs)
	    run.regenerateOutputFiles(cl.getCheckpointId());
	}

#ifdef _DISABLE_INIT_LNL_CHECK
      tout << std::endl << "DEVEL: I will not check, if resuming the chains yields the exact same\n" 
	   << "likelihood (that we had when writing the checkpoint). For sequential\n"
	   << "runs, this should be fine, but if you built with MPI, you the same\n"
	   << "seed may not result in the exact same run" << std::endl;       
#endif
    }

  tout << endl << "Will run " << runParams.getNumRunConv() << " runs in total with "<<  runParams.getNumCoupledChains() << " coupled chains" << endl << endl;   

  // determine and print initial state for each chain  
  tout << std::endl << "initial state: " << endl; 
  tout << "================================================================" << endl; 
  for(auto &run: runs)
    {
      for(auto &chain :  run.getChains())
	{
	  chain.resume(true, 
#ifdef _DISABLE_INIT_LNL_CHECK
		       false
#else 
		       true 
#endif
		       ); 
	  chain.suspend(false);
	  
	  tout << chain; 
	  tout << "\tRNG(" << chain.getChainRand() << ")"<< endl; 	  
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
	  ss <<  PROGRAM_NAME << "_topologies." << cl.getRunid() << "." << i; 
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
	    tout << std::endl  << "ASDSF for trees " << asdsf.getStart() << "-" << asdsf.getEnd() << ": " <<setprecision(2) << asdsfVal * 100 << "%" << std::endl << std::endl; 

	  return asdsfVal < runParams.getAsdsfConvergence(); 

	}
      else 
	return false; 
      
    }
  else 
    return true; 
}


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
				      vector<unique_ptr<AbstractParameter> > &variableResult, 
				      shared_ptr<LikelihoodEvaluator> eval)
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


void SampleMaster::printDuringRun(nat gen)
{
  tout << "[" << gen << "," << CLOCK::duration_cast<CLOCK::duration<double> > (CLOCK::system_clock::now()- timeIncrement   ).count()     << "s]\t"; 
  timeIncrement = CLOCK::system_clock::now(); 

  bool isFirst = true; 
  for(auto &run : runs)
    {
      if(isFirst)
	isFirst = false; 
      else 
	tout << " === " ; 

      auto &chains = run.getChains();
      std::vector<const Chain*> chainPtrs; 
      for(auto &c : chains)
	chainPtrs.push_back(&c);       
      sort(chainPtrs.begin(), chainPtrs.end(), [] (const Chain* a , const Chain* b ) {return a->getCouplingId() < b->getCouplingId(); }); 

      bool isFirstL = true;  
      for(auto &c: chainPtrs)
	{
	  if(isFirstL)
	    isFirstL = false; 
	  else 
	    tout << " "; 
	  tout << c->getLikelihood(); 	  
	}
    }  
  tout << std::endl; 
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

  bool hasConverged = false;   
  int curGen = runs[0].getChains()[0].getGeneration(); 
  
  int lastPrint = (curGen / runParams.getPrintFreq() )  * runParams.getPrintFreq() ; 
  int lastDiag =  (curGen / runParams.getDiagFreq() ) * runParams.getDiagFreq(); 
  int lastChkpnt = ( curGen / runParams.getChkpntFreq() )* runParams.getChkpntFreq() ; 
  
  printDuringRun(curGen); 

  while(curGen < runParams.getNumGen() || not hasConverged)   
    { 
      int nextPrint =  lastPrint + runParams.getPrintFreq(); 
      int nextDiag = lastDiag + runParams.getDiagFreq(); 
      int nextChkpnt = lastChkpnt + runParams.getChkpntFreq(); 
      int toExecute = min(nextChkpnt, 
			  min(nextPrint, nextDiag)) -  curGen; 

      for(auto &run : runs)
	run.executePart(toExecute, pl );

      curGen += toExecute; 


      if(curGen % runParams.getDiagFreq() == 0 )
	{
	  hasConverged = convergenceDiagnostic(); 
#ifdef ENABLE_PRSF      
	  printPRSF(run_id);
#endif
	  lastDiag = curGen; 
	}
	
      if(curGen % runParams.getPrintFreq() == 0 )
	{
	  printDuringRun(curGen); 
	  lastPrint = curGen; 
	}

      if(curGen % runParams.getChkpntFreq() == 0 )
	{
	  tout << "writing checkpoint" << std::endl; 
	  writeCheckpointMaster();
	  lastChkpnt = curGen; 
	} 
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


void SampleMaster::readFromCheckpoint( std::ifstream &in ) 
{
  for(auto &run : runs)
    run.readFromCheckpoint(in); 

  long start = 0; 
  in >> start; 
  CLOCK::duration<long> durationSinceStart{start};   
  initTime = CLOCK::system_clock::now() -  durationSinceStart; 
}
 
void SampleMaster::writeToCheckpoint( std::ofstream &out) 
{  
  for(auto & run : runs)    
    run.writeToCheckpoint(out);

  long duration = CLOCK::duration_cast<CLOCK::duration<long> >(CLOCK::system_clock::now() - initTime).count();
  cWrite(out, duration); 
}


void SampleMaster::writeCheckpointMaster()
{  
  if(pl.isMasterReporter() )
    {
      stringstream ss; 
      ss <<  PROGRAM_NAME << "_newCheckpoint." << cl.getRunid()  ; 
      std::string newName = ss.str();
      ofstream chkpnt; 
      Checkpointable::getOfstream(newName, chkpnt); 
      writeToCheckpoint(chkpnt);
      chkpnt.close(); 
	  
      ss.str("");
      ss << PROGRAM_NAME << "_checkpoint." << cl.getRunid(); 
      std::string curName = ss.str();
      if( std::ifstream(curName) )
	{
	  ss.str("");
	  ss << PROGRAM_NAME << "_prevCheckpointBackup." << cl.getRunid(); 
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
      std::ofstream nullstream("/dev/zero"); 
      writeToCheckpoint(nullstream); 
      nullstream.close(); 
    }
} 


#include "IntegrationModuleImpl.hpp"
