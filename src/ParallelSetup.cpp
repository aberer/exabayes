#include <iostream> 
#include <sstream>

#include "ParallelSetup.hpp"
#include "GlobalVariables.hpp"
#include "SwapMatrix.hpp"
#include "SampleMaster.hpp"


ParallelSetup::ParallelSetup(int argc, char **argv)
  : runsParallel(1)
  , chainsParallel(1)
  , globalRank(0)
  , globalSize(1)
{  
#if HAVE_PLL == 0 
  MPI_Init(&argc, &argv);
#endif
}




void ParallelSetup::printLoadBalance(const TreeAln& traln ) const 
{
#if HAVE_PLL == 0 
  auto &&ss = std::stringstream{};

  ss << SOME_FIXED_PRECISION; 

  if(isGlobalMaster())
    ss << "load distribution (rank,#pattern):" << std::endl; 
  
  double sumA = 0;   

  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      auto& partition = traln.getPartition(i); 
      sumA += partition.width; 
    }
  
  if(isGlobalMaster())
    ss << "[ " << globalRank << " ] " ; 
  ss  << globalRank<< "\t" << sumA << std::endl; 

  if(sumA == 0)
    {
      ss << std::endl << "WARNING: there is a process that does not evaluate any portion of the\n"
	 << "data. Consider enabling or disabling the '-Q' option"  << std::endl; 
    }
  
  blockingPrint(MPI::COMM_WORLD, ss.str() );
  MPI_Barrier(MPI::COMM_WORLD); 
#endif
}


// TODO asserts with generations 


void ParallelSetup::synchronizeChainsAtMaster( vector<CoupledChains>& runs, CommFlag commFlags) const
{  
  std::stringstream ss; 
#if 0   
  ss << "Communicating " ; 
  if(commFlags & CommFlag::Swap)
    ss << " SWAP, " ; 
  if(commFlags & CommFlag::PrintStat)
    ss << " PRINTSTAT, "; 
  if(commFlags & CommFlag::Tree)
    ss << " TREE, "; 
  if(commFlags & CommFlag::Proposals)
    ss << " PROPOSALS, "; 
  ss << std::endl; 
  blockingPrint(chainLeaderComm, ss.str()); 
  ss.str(""); 
#endif

  std::vector<char> serialized; 

#if HAVE_PLL == 0  
  // serializing all chains that are assigned to me  
  int lengthOfChainStream = 0; 
  bool isFirst = true ; 
  for(auto &run : runs)
    {
      for(nat i = 0; i < run.getChains().size() ;++i)
	{
	  const Chain& chain = run.getChains()[i]; 
	  if(isMyChain(run.getRunid(), i))
	    {
	      std::string chainSer = chain.serializeConditionally( commFlags ); 

	      if(isFirst)
		lengthOfChainStream = chainSer.size(); 
	      isFirst = false; 
	      assert(nat(lengthOfChainStream) == chainSer.size()); 

	      serialized.insert(serialized.end(), chainSer.begin(), chainSer.end()); 
	    }
	}
    }
  assert(not isFirst); 

  int lengthOfSwap = 0;   
  if(commFlags & CommFlag::Swap)
    {
      isFirst = true; 
      for(auto &run : runs)
	{
	  if(isMyRun(run.getRunid()))
	    {
	      std::stringstream part; 
	      run.getSwapInfo().serialize(part); 
	      auto asString = part.str();
	      
	      if(isFirst )
		lengthOfSwap = asString.size(); 

	      isFirst = false; 
	      assert(lengthOfSwap == int(asString.size())); 
	      serialized.insert(serialized.end(), asString.begin(), asString.end()); 
	    }
	}
      assert(not isFirst); 
    }

#ifdef EFFICIENT
  // we communicate the swap matrix multiple times  => gatherv necessary 
  assert(0); 
#endif

  std::string allMyChains (serialized.begin(), serialized.end());   
  assert(allMyChains.size( )== serialized.size()); 
  nat lengthPerProcess = allMyChains.size(); 

  nat totalLength = chainLeaderComm.Get_size() * lengthPerProcess; 
  std::vector<char> allChainsSerialized;
  allChainsSerialized.resize(totalLength); 

  // communicate 
  chainLeaderComm.Gather( allMyChains.c_str(), lengthPerProcess, MPI::CHAR, &(allChainsSerialized[0]), lengthPerProcess, MPI::CHAR, 0 );

  // sorry for all that nesting...
  auto streamIter = allChainsSerialized.begin(); 
  if(chainLeaderComm.Get_rank() == 0)
    {
      for(int procIdx  = 0; procIdx < chainLeaderComm.Get_size(); ++procIdx)
	{
	  nat rankOfChainLeader = procIdx * getProcessesPerChainBatch(); 

	  for(auto &run : runs)
	    {	      
	      for(nat i = 0; i < run.getChains().size(); ++i)
		{		  
		  if(isMyChain(rankOfChainLeader, run.getRunid(), i))
		    {		      
		      auto &chain = run.getChains().at(i);
		      std::string str(streamIter, streamIter + lengthOfChainStream ); 
		      streamIter += lengthOfChainStream; 
		      chain.deserializeConditionally(str, commFlags); 
		    }
		}
	    }
	  
	  if(commFlags & CommFlag::Swap)
	    {
#ifdef UNSURE
	      // this is not perfect: we only get the global acc/rej
	      // values but it is really not worth it, spending much
	      // more thought on this...
	      assert(0); 
#endif

	      for(auto &run : runs)
		{
		  if(isMyRun(rankOfChainLeader, run.getRunid()))
		    { 
		      ss.str() = std::string(streamIter, streamIter + lengthOfSwap); 
		      streamIter += lengthOfSwap; 

		       auto sw = SwapMatrix(run.getChains().size()); 
		      sw.deserialize(ss); 
		      // randCtr_t v = {0,0} ; 		      
		      // Randomness rand(v);  
		      // rand.readFromCheckpoint(ss);
		      // std::cout << "READ Randomness " << rand << std::endl; 
		      // run.setRandomness(rand); 

		      // HACK: with more efficient communication, we
		      // would need to do this only once
		      if(isRunLeader(procIdx))
			run.addToSwapMatrix(sw); 		      
		    }
		}
	    }
	}
    }
#endif  

  // std::cout << << "made it!!"
  // ss.str("") ;
  // ss << "made it" << std::endl; 
  // blockingPrint(chainLeaderComm, ss.str()); 
}

void ParallelSetup::finalize()
{
#if HAVE_PLL == 0
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
}



void ParallelSetup::genericExit(int code)
{
#if HAVE_PLL == 0
  // MPI_(MPI_COMM_WORLD, code); 
  MPI_Finalize(); 
  exit(code); 
#else 
  exit(code);   
#endif
}



std::ostream& operator<<(std::ostream& out, const ParallelSetup &pl)
{
  out << "[GLOBAL:" << pl.globalRank << "/" << pl.globalSize
  << ",RUN:" <<   pl.getRankInRunBatch() <<  "/" << pl.getProcessesPerRunBatch()
  << ",CHAIN:" << pl.getRankInChainBatch() << "/" << pl.getProcessesPerChainBatch()
  <<  "]"; 
  return out; 
}


bool 
ParallelSetup::isMyChain(nat gRank, nat runid, nat chainIndex) const 
{
  return  isMyRun(gRank, runid) && (chainIndex % chainsParallel == getChainBatch(gRank)); 
}


bool 
ParallelSetup::isMyChain(nat runid, nat chainIndex) const 
{
  return isMyChain(globalRank, runid, chainIndex); 
}


bool 
ParallelSetup::isMyRun(nat runid) const 
{
  return isMyRun(globalRank, runid); 
}

bool 
ParallelSetup::isMyRun(nat gRank, nat runid) const 
{
  return runid % getRunsParallel() == getRunBatch(gRank);
}



std::string
ParallelSetup::sendRecvChain(const CoupledChains& run, nat myIndex, nat otherChainIndex, std::string myChainSer, CommFlag flags )  
{
  // trivial, if we own the chain 
  if(isMyChain(run.getRunid(), otherChainIndex))
    {
      // std::cout << "both mine" << std::endl; 
      return run.getChains()[otherChainIndex].serializeConditionally(flags); 
    }
  else 
    {      
      std::vector<char> otherChainSer; 
      otherChainSer.resize(myChainSer.size());

#if HAVE_PLL == 0
      nat otherChainBatch = otherChainIndex % chainsParallel; 
      
      if(myIndex < otherChainIndex) 
	{
	  // send first 	  
	  commToChains[otherChainBatch].Bcast(&(myChainSer[0]), myChainSer.size(), MPI::CHAR, isChainLeader() ? MPI_ROOT : MPI_PROC_NULL );
	  
	  // now receive
	  commToChains[otherChainBatch].Bcast(&(otherChainSer[0]), myChainSer.size(), MPI::CHAR, 0); 
	}
      else 
	{
	  // receive first
	  commToChains[otherChainBatch].Bcast(&(otherChainSer[0]), myChainSer.size(), MPI::CHAR, 0); 

	  // now send 	  
	  commToChains[otherChainBatch].Bcast(&(myChainSer[0]), myChainSer.size(), MPI::CHAR, isChainLeader() ? MPI_ROOT : MPI_PROC_NULL );
	}
#else 
      assert(0);
#endif

      std::string result(otherChainSer.begin(), otherChainSer.end()); 
      return result; 
    }

}


#if HAVE_PLL == 0 
void
ParallelSetup::blockingPrint(const MPI::Comm &comm,std::string ss )
{
  int myRank = comm.Get_rank(); 
  int size = comm.Get_size(); 

  if(myRank != 0)
    {
      int res = 0; 

      
      comm.Recv(&res, 1, MPI::INT, myRank-1, 0);       
      assert(res == 1 ); 
    }

  std::cout << "[ " << myRank << " ] " << ss << std::endl; 
  
  if(myRank != size - 1) 
    {
      int res = 1; 
      comm.Send(&res, 1 , MPI::INT, myRank + 1, 0);       
    }  
}


void ParallelSetup::initializeExaml(const CommandLine &cl)
{   
  globalRank = MPI::COMM_WORLD.Get_rank(); 
  globalSize = MPI::COMM_WORLD.Get_size(); 

  runsParallel = cl.getNumRunParallel();
  chainsParallel = cl.getNumChainsParallel();

  if(globalSize <  runsParallel)
    {
      if(globalRank == 0)
	{
	  std::cout << "You requested to run "  << runsParallel << " in parallel, however there are only " 
		    << globalSize << " processes (we need at least 1 process per run, see command line option -R )"  << std::endl; 
	}
      MPI::COMM_WORLD.Abort(1); 
    }  

  if(globalSize % (runsParallel * chainsParallel) != 0) 
    {
      std::cerr << "You called " << PROGRAM_NAME << " with " 
  		<< globalSize << " processes and instructed " << PROGRAM_NAME
  		<< " to run " << runsParallel << " independent runs in parallel.\n"
  		<< "For efficiency, please make sure that the number of procesess\n"
  		<< "is a multiple of the product of runs to be executed in\n" 
  		<< "parallel and the number of chains to be executed in"
  		<< "parallel. Aborting." << std::endl; 
      exit(0);
    }

  runComm = MPI::COMM_WORLD.Split( getRunBatch(), getRankInRunBatch() ); 
  chainComm = runComm.Split(getChainBatch(), getRankInChainBatch() ); 

  // setup intecomms between chains 
  for(nat i = 0; i < chainsParallel; ++i )
    {
      nat myChainBatch =  getChainBatch() ; 
      if(i == getChainBatch())
	continue; 
      
      int tag = std::min(i, myChainBatch) * chainsParallel + std::max(i, myChainBatch); 
      int rLeader = getRunBatch() * getProcessesPerRunBatch() + i * getProcessesPerChainBatch() ; 

      std::stringstream ss; 
      ss << "INTERCOMM CREATE: contacting " << rLeader << " with tag " << tag << std::endl; 
      blockingPrint(runComm, ss.str()); 
      
      commToChains[i] = chainComm.Create_intercomm(0, MPI::COMM_WORLD, rLeader, tag); 
    }

  std::stringstream ss; 
  if(isChainLeader())
    ss << "is Chain Leader" << std::endl; 
  blockingPrint(MPI::COMM_WORLD, ss.str()); 


  // do it again for the axml communicator 
  MPI_Comm tmpComm;   
  MPI_Comm_split(MPI_COMM_WORLD, getRunBatch(), getRankInRunBatch(), &tmpComm); 
  MPI_Comm_split(tmpComm, getChainBatch(), getRankInChainBatch(), &comm);   
  MPI_Comm_free(&tmpComm); 
  MPI_Comm_rank(comm, &processID); 
  MPI_Comm_size(comm, &processes); 

  // initialize the comm connection all chain leaders 
  chainLeaderComm = MPI::COMM_WORLD.Split(getRankInChainBatch(), getRunBatch() * getProcessesPerRunBatch() + getChainBatch()  ); 

  // print your configuration     
  ss.str(""); 
  ss << *this << std::endl; 
  blockingPrint(MPI::COMM_WORLD,  ss.str() ); 
}
#endif
