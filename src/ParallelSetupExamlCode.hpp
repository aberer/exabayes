#include <iostream> 
#include <sstream>
#include "GlobalVariables.hpp"
#include "SwapMatrix.hpp"
#include "SampleMaster.hpp"
#include "IncompleteMesh.hpp"

#include "common.h"


ParallelSetup::ParallelSetup()
  : _runsParallel(1)
  , _chainsParallel(1)
{  
  _mesh.initialize(1,1,1);
}


ParallelSetup::ParallelSetup(ParallelSetup&& rhs)
  : _runsParallel(std::move(rhs._runsParallel))
  , _chainsParallel(std::move(rhs._chainsParallel))
  , _globalComm{  std::move(rhs._globalComm)}
  , _runComm{ std::move(rhs._runComm)}
  , _chainComm { std::move(rhs._chainComm)}
  , _commToChains{ std::move(rhs._commToChains)}
  , _mesh(std::move(rhs._mesh))
  , _chainLeaderComm(std::move(rhs._chainLeaderComm))
{
}


ParallelSetup& ParallelSetup::operator=(ParallelSetup rhs)
{
  swap(*this, rhs); 
  return *this; 
}


void swap(ParallelSetup &lhs, ParallelSetup &rhs)
{
  using std::swap; 
  swap(lhs._runsParallel, rhs._runsParallel); 
  swap(lhs._chainsParallel, rhs._chainsParallel); 
  swap(lhs._runComm, rhs._runComm); 
  swap(lhs._chainComm, rhs._chainComm);; 
  swap(lhs._commToChains, rhs._commToChains);
  swap(lhs._mesh, rhs._mesh); 
  swap(lhs._globalComm, rhs._globalComm); 
  swap(lhs._chainLeaderComm, rhs._chainLeaderComm); 
}


void ParallelSetup::initialize(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
}


void ParallelSetup::printLoadBalance(const TreeAln& traln, nat numRun, nat numCoupled ) const 
{
  auto &&ss = std::stringstream{};

  ss << SOME_FIXED_PRECISION; 

  if(isGlobalMaster())
    ss << "load distribution (rank,#pattern,chainsPerRun):"  << std::endl;
  
  double sumA = 0;   

  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      auto& partition = traln.getPartition(i); 
      sumA += partition.width; 
    }
  
  if(isGlobalMaster())
    ss << "[ " << _globalComm.getRank() << " ] " ; 
  ss  << sumA << "\t"  ; 

  bool isFirst = true; 
  for(nat i = 0; i < numRun; ++i )
    {
      if(not isMyRun(i))
	continue; 
      for(nat j = 0; j < numCoupled; ++j)
	{
	  if(not isMyChain(i,j))
	    continue; 
	  ss << (not isFirst ? "," : "") << "("<< i << "," << j << ")" ; 
	  isFirst = false; 
	}
    }

  ss << std::endl; 
  
  if(sumA == 0)
    {
      ss << std::endl << "WARNING: there is a process that does not evaluate any portion of the\n"
	 << "data. Consider enabling or disabling the '-Q' option"  << std::endl; 
    }

  auto tmp = Communicator{}; 
  blockingPrint(tmp , ss.str() );
}




auto ParallelSetup::serializeAllChains(std::vector<CoupledChains> &runs, CommFlag commFlags) const
  -> std::tuple<nat,std::vector<char>> 
{
  auto serialized = std::vector<char>{}; 
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
	      auto &&ss = std::ostringstream{}; 
	      chain.serializeConditionally( ss, commFlags ); 
	      auto chainSer =  ss.str();

	      if(isFirst)
		lengthOfChainStream = chainSer.size(); 
	      isFirst = false; 
	      assert(nat(lengthOfChainStream) == chainSer.size()); 

	      serialized.insert(end(serialized), begin(chainSer), end(chainSer)); 
	    }
	}
    }

  // TODO 
  assert(not isFirst); 

  return std::make_tuple(lengthOfChainStream,serialized); 
}



auto ParallelSetup::serializeSwapMatrices(std::vector<CoupledChains>& runs, CommFlag &commFlags ) const 
  -> std::tuple<nat, std::vector<char>> 
{
  auto serialized = std::vector<char>{};
  int lengthOfSwap = 0;   
  if( ( commFlags & CommFlag::SWAP ) != CommFlag::NOTHING )
    {
      bool isFirst = true; 
      for(auto &run : runs)
	{
	  if(isMyRun(run.getRunid()))
	    {
	      auto && part = std::stringstream{}; 
	      run.getSwapInfo().serialize(part); 

	      auto asString = part.str();
	      
	      if(isFirst)
		lengthOfSwap = asString.size();
	      
	      isFirst = false; 
	      assert(lengthOfSwap == int(asString.size())); 
	      serialized.insert(end(serialized), begin(asString), end(asString)); 
	    }
	}
      assert(not isFirst); 
    }

  return std::make_tuple(lengthOfSwap, serialized);
}




// TODO i guess we communicate the swap matrix mutliple times, this could be saved
void ParallelSetup::synchronizeChainsAtMaster( std::vector<CoupledChains>& runs, CommFlag commFlags) const
{

#if 0   
  ss << "Communicating " ; 
  if( ( commFlags & CommFlag::Swap ) != CommFlag::NOTHING)
    ss << " SWAP, " ; 
  if( ( commFlags & CommFlag::PrintStat )  != CommFlag::NOTHING)
    ss << " PRINTSTAT, "; 
  if( ( commFlags & CommFlag::Tree ) != CommFlag::NOTHING)
    ss << " TREE, "; 
  if( ( commFlags & CommFlag::Proposals ) != CommFlag::NOTHING)
    ss << " PROPOSALS, "; 
  ss << std::endl; 
  blockingPrint(chainLeaderComm, ss.str()); 
  ss.str(""); 
#endif
  
  auto myData = std::vector<char>{}; 
  auto lengthOfChainStream = 0; 
  std::tie(lengthOfChainStream, myData) = serializeAllChains(runs, commFlags );

  nat lengthOfSwap = 0; 
  auto swapSer = std::vector<char>{}; 
  std::tie(lengthOfSwap, swapSer) = serializeSwapMatrices(runs, commFlags);
  myData.insert(end(myData), begin(swapSer), end(swapSer)); 

  // communicate 
  auto allChainsSerialized = _chainLeaderComm.gather(myData);
   

  // the master process unpacks the received information 
  auto streamIter = begin(allChainsSerialized); 
  if( isGlobalMaster() )
    {
      auto ranksofLeaders = std::vector<nat> {}; 
      for(int i = 0; i < _globalComm.getSize(); ++i)
	{
	  auto coords =  _mesh.getCoordinates(i);
	  if(coords[2] == 0) 
	    ranksofLeaders.push_back(i); 
	}

      assert(ranksofLeaders.size() == nat(_chainLeaderComm.getSize())); 

      // accumulate swap matrices 
      auto sws = std::vector<SwapMatrix>{}; 
      for(nat i = 0; i < runs.size() ;++i)
	sws.emplace_back(runs[0].getChains().size());

      for(auto rankOfChainLeader : ranksofLeaders)
	{
	  for(auto &run : runs)
	    {
	      for(nat i = 0; i < run.getChains().size(); ++i)
		{		  
		  if(isMyChain(rankOfChainLeader, run.getRunid(), i))
		    {		      
		      auto &chain = run.getChains().at(i);
		      auto &&iss =  std::istringstream{std::string(streamIter, streamIter + lengthOfChainStream )}; 
		      streamIter += lengthOfChainStream; 
		      chain.deserializeConditionally(iss, commFlags); 
		    }
		}
	    }
	  
	  if( ( commFlags & CommFlag::SWAP ) != CommFlag::NOTHING )
	    {
	      for(auto &run : runs)
		{
		  if(isMyRun(rankOfChainLeader, run.getRunid()))
		    { 
		      auto&& iss = std::istringstream{std::string{streamIter, streamIter + lengthOfSwap}};
		      streamIter += lengthOfSwap; 
		      auto sw = SwapMatrix(run.getChains().size()); 
		      sw.deserialize(iss); 

		      // HACK: with more efficient communication, we
		      // would need to do this only once

		      if(isRunLeader(rankOfChainLeader))
			sws[run.getRunid()] += sw ; 
		    }
		}
	    }
	}
      
      if( ( commFlags & CommFlag::SWAP )  != CommFlag::NOTHING)
	{
	  // set the accumulated swap matrices 
	  for(auto &run : runs)
	    run.setSwapInfo(sws.at(run.getRunid()));
	}
    }
  else
    {
      // for reproducibility. Otherwise runs are not reproducible any
      // more with chain-level parallelism
      if(  ( commFlags & CommFlag::TREE ) != CommFlag::NOTHING )
      	{
      	  for(auto &run : runs)
      	    {
      	      nat ctr = 0; 
      	      for(auto &chain : run.getChains())
      		{
      		  if(isMyChain(run.getRunid(), ctr))
      		    {
      		      auto &&oss = std::ostringstream{};
      		      chain.serializeConditionally(oss, commFlags );
      		      auto &&iss = std::istringstream{oss.str()}; 
      		      chain.deserializeConditionally(iss,  commFlags);
      		    }
      		  ++ctr;
      		}
      	    }
      	}
    }
}


void ParallelSetup::finalize()
{
  MPI_Finalize();
}



void ParallelSetup::genericExit(int code)
{
  if( code == 0 )
    MPI_Finalize(); 
  else 
    {
      MPI_Abort(MPI_COMM_WORLD, code);
    }
  exit(code); 
}


std::ostream& operator<<(std::ostream& out, const ParallelSetup &pl)
{
  auto coords = pl._mesh.getCoordinates(pl.getGlobalRank()); 
  out << "[" << coords[0] << ","<< coords[1] << "," << coords[2]  << "]"; 
  return out; 
}


bool 
ParallelSetup::isMyChain(nat gRank, nat runid, nat chainIndex) const 
{
  auto hisCoords= _mesh.getCoordinates(gRank);
  return isMyRun(gRank, runid ) &&  ( chainIndex % getChainsParallel()) == hisCoords[1]; 
}


bool 
ParallelSetup::isMyChain(nat runid, nat chainIndex) const 
{
  return isMyChain(getGlobalRank(), runid, chainIndex); 
}


bool 
ParallelSetup::isMyRun(nat runid) const 
{
  return isMyRun(getGlobalRank(), runid);
}

bool 
ParallelSetup::isMyRun(nat gRank, nat runid) const 
{
  auto hisCoords = _mesh.getCoordinates(gRank);
  return (runid % getRunsParallel() )  == hisCoords[0]; 
}


std::string
ParallelSetup::sendRecvChain(const CoupledChains& run, nat myIndex, nat otherChainIndex, std::string myChainSer, CommFlag flags )  
{
  nat runid = run.getRunid();
  assert(isMyChain(runid, myIndex) && not isMyChain(runid, otherChainIndex)); 

  auto otherChainSer = std::vector<char>{} ; 
  otherChainSer.resize(myChainSer.size());

  nat otherChainBatch = otherChainIndex % _chainsParallel; 

  if(myIndex < otherChainIndex) 
    {
      // send first 
      MPI_Bcast((void*)myChainSer.data(), myChainSer.size(), MPI_CHAR, isChainLeader() ? MPI_ROOT : MPI_PROC_NULL , _commToChains[otherChainBatch].getHandle()); 

      // now receive
      MPI_Bcast((void*)otherChainSer.data(), myChainSer.size(), MPI_CHAR, 0, _commToChains[otherChainBatch].getHandle()); 
    }
  else 
    {
      // receive first
      MPI_Bcast((void*)otherChainSer.data(), myChainSer.size(), MPI_CHAR, 0, _commToChains[otherChainBatch].getHandle());

      // now send 	  
      MPI_Bcast((void*)myChainSer.data(), myChainSer.size(), MPI_CHAR, isChainLeader() ? MPI_ROOT : MPI_PROC_NULL, _commToChains[otherChainBatch].getHandle() ); 
    }
  
  auto result = std::string (otherChainSer.begin(), otherChainSer.end()); 
  return result; 
}


void ParallelSetup::blockingPrint( Communicator &comm,std::string ss )
{
  int myRank = comm.getRank();
  int size = comm.getSize(); 

  if(myRank != 0)
    {
      int res = comm.receive<int>(myRank-1,0);
      assert(res == 1 ); 
    }

  std::cout << "[ " << myRank << " ] " << ss << std::endl; 
  std::cout.flush();
  
  if(myRank != size - 1) 
    comm.send<int>(1, myRank + 1 , 0 ); 
  else 
    comm.send<int>(1,0,0 ); 
    
  if( myRank == 0  )
    {
      auto res =  comm.receive<int>(size-1,0);
      assert(res == 1); 
    }
}


void ParallelSetup::initializeExaml(const CommandLine &cl )
{
  _globalComm = Communicator{}; 
  assert(_globalComm.isValid()) ; 
  assert(_globalComm.getSize() > 0 ); 
  
  _runsParallel = cl.getNumRunParallel();
  _chainsParallel = cl.getNumChainsParallel();

  if( _globalComm.getSize() < int(_runsParallel * _chainsParallel)  && _globalComm.getRank() == 0 )
    {
      std::cout << "\nError: You attempted to run " << _runsParallel << " runs in parallel for which " <<  _chainsParallel <<  " coupled chains should be run in parallel. This requires at least " << (_runsParallel * _chainsParallel) << " processes (in the best case a multiple of this number). " << PROGRAM_NAME << " was started with " << _globalComm.getSize() << " processes though.\n" <<  std::endl; 
      genericExit(-1);
    }

  if( _globalComm.getSize() % (_runsParallel * _chainsParallel) != 0  && _globalComm.getRank() == 0 ) 
    {
      std::cout << "\tWARNING! The number of processes and number chains run in parallel\n"
		<< "\tare not a multiple of the overall nmuber of processes. Parts of the\n"
		<< "\tcode may be executed less efficiently.\n" ; 
    }

  _mesh.initialize(_globalComm.getSize(), _runsParallel, _chainsParallel);

  auto coords = _mesh.getCoordinates(_globalComm.getRank());

  // create the run comm 
  _runComm = Communicator().split(coords[0], coords[1] * coords[2]);
  assert(_runComm.isValid()); 
  
  // create chain comm 
  _chainComm = _runComm.split(coords[1], coords[2]);
  assert(_chainComm.isValid()); 
  
  // create chain intercomms  
  for(nat i = 0; i < _chainsParallel; ++i)
    {
      nat myRunBatch = coords[0]; 
      nat myChainBatch = coords[1] ; 
      if(i == myChainBatch)
	continue; 
      
      int tag = std::min(i, myChainBatch) * _chainsParallel + std::max(i ,myChainBatch); 
      auto tmp = std::array<nat,3>{{myRunBatch, i , 0 }}; 
      int rLeader = _mesh.getRankFromCoordinates( tmp );

      // auto &&ss = std::stringstream{}; 
      // auto myCoords = getMyCoordinates();
      // ss << "[" <<  getGlobalRank()  << "] INTERCOMM CREATE: contacting " << rLeader << " with tag " << tag << " my batch: " << myChainBatch << "," << myRunBatch << " and coords " << myCoords[0] << "," << myCoords[1] << "," << myCoords[2] << std::endl; 
      // std::cout << ss.str() << std::endl; 
      // blockingPrint(_runComm, ss.str()); 

      MPI_Intercomm_create(_chainComm.getHandle(), 0, MPI_COMM_WORLD, rLeader, tag, _commToChains[i].getPtr());

      assert(_commToChains[i].isValid()); 
    } 

  // auto && ss = std::stringstream{} ; 
  // if(not cl.isQuiet() )
  //   {
  //     if(coords[2] == 0)
  // 	ss << "is Chain Leader" << std::endl; 
  //     blockingPrint(_globalComm, ss.str()); 
  //   }
  
  // install the ExaML communicator 
  MPI_Comm_dup(_chainComm.getHandle(), &comm); 
  assert(comm != MPI_COMM_NULL); 
  MPI_Comm_rank(comm, &processID);
  MPI_Comm_size(comm, &processes); 

  // this can be afforded once, but it is important, that everything
  // is mapped accurately
  nat myRank = 0; 
  for(int i = 0; i < _globalComm.getSize() ; ++i)
    {
      if(i == _globalComm.getRank())
	break; 
      auto hisCoords = _mesh.getCoordinates(i);
      if( hisCoords[2] == coords[2])
	++myRank;
    }

  // create chain leader comm 
  _chainLeaderComm = _globalComm.split(coords[2], myRank); 

  // if(not cl.isQuiet() )
  //   {
  //     ss.str(""); 
  //     ss << *this << std::endl; 
  //     blockingPrint(_globalComm,  ss.str() ); 
  //   }
}

nat ParallelSetup::getRankInData() const 
{
  auto coords = _mesh.getCoordinates(getGlobalRank());

  auto result = coords[2]; 
  assert(int(result) == _chainComm.getRank()); 
  return result; 
}


nat ParallelSetup::getSizeInData() const 
{
  return _chainComm.getSize(); 
}


bool ParallelSetup::isRunLeader(nat gRank) const 
{  
  auto coords = _mesh.getCoordinates(getGlobalRank());
  return coords[1] == 0 && coords[2] == 0; 
}  


bool ParallelSetup::isChainLeader() const
{
  auto coords = _mesh.getCoordinates(getGlobalRank()); 
  return coords[2] == 0; 
}


bool ParallelSetup::isGlobalMaster() const 
{
  auto coords = _mesh.getCoordinates(getGlobalRank()); 
  return coords[0] == 0 
    && coords[1] == 0
    && coords[2] == 0; 
}


bool ParallelSetup::globalBroadcast( bool val,int root ) const
{
  return _globalComm.broadcast(val, root);
}



nat ParallelSetup::chainIdToLeaderRank(nat runId, nat chainId) const  
{
  nat pBatch = runId % _runsParallel; 
  nat cBatch = chainId % _chainsParallel; 
  
  auto coords = std::array<nat,3> {{ pBatch, cBatch, 0 }}; 
  return _mesh.getRankFromCoordinates(coords);
}


void ParallelSetup::globalBarrier()  
{
  _globalComm.waitAtBarrier(); 
}


int ParallelSetup::getRankInRun( std::array<nat,3> coords) const 
{
  auto localRootCoords = std::array<nat,3> { { coords[0], 0, 0  }}; 
  return _mesh.getRankFromCoordinates(coords)  - _mesh.getRankFromCoordinates(localRootCoords); 
} 
