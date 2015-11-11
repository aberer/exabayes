#ifndef _PARALLEL_IMPL_EXAML_HPP
#define _PARALLEL_IMPL_EXAML_HPP

#include <unordered_map>
#include "CoupledChains.hpp"
#include "config/CommandLine.hpp"
#include "IncompleteMesh.hpp"
#include "FlagType.hpp" 
#include "Communicator.hpp"


#include <mpi.h>
extern MPI_Comm comm; 		// the examl communicator 
extern int processID; 		// examl rank 
extern int processes; 		// examl comm size 

class ParallelSetup
{
public: 
  ParallelSetup();
  ParallelSetup(ParallelSetup &&rhs); 
  ParallelSetup(const ParallelSetup& rhs)  = delete ; 
  ParallelSetup& operator=( ParallelSetup rhs) ; 
  friend void swap(ParallelSetup &lhs, ParallelSetup &rhs); 

  static void initialize(int argc, char **argv); 
  
  /** 
      @brief finalize the parallel environment 
   */ 
  static void finalize(); 
  /** 
      @brief gets the number of runs executed in parallel  
   */ 
  nat getRunsParallel() const {return _runsParallel ; }
  nat getGlobalRank() const { return _globalComm.getRank(); }
  nat getGlobalSize() const { return _globalComm.getSize(); }
  /** 
      @brief gets the number of (coupled) chains executed in paralel (in addition to run-level parallelism)
   */ 
  nat getChainsParallel() const {return _chainsParallel; }

  nat getRankInData() const ; 
  nat getSizeInData() const ; 
  /** 
      @brief indicates whether the process should conduct output operatons (for his run) 
   */ 
  bool isRunLeader() const {  return isRunLeader( getGlobalRank() ) ; }  
  bool isRunLeader(nat gRank) const ; 
  /**
     @brief indicates whether a process is the leader in a batch of chains 
   */ 
  bool isChainLeader() const; 
  /** 
      @brief indicates whether the process is the process to conduct global logging output
   */ 
  bool isGlobalMaster() const ; 
  /** 
      @brief sends a serialized chain representation to the other involved process  
      @param run relevant run  
      @param otherChainIndex the index of the chain that is involved in the swapping attempt   
      @return serialized representation  of the chain we swap with 
   */ 
  std::string
  sendRecvChain(const CoupledChains& run, nat myIndex, nat otherChainIndex  , std::string myChainSer, CommFlag flags )  ; 

  nat chainIdToLeaderRank(nat runId, nat chainId) const ; 

  /** 
      @brief synchronize all core chain information at the master node
      @notice this only concerns the global masterprocess in a
      parallel environtment.  
  */ 
  void synchronizeChainsAtMaster( std::vector<CoupledChains>& runs, CommFlag commFlags) const;  

  /** 
      @brief initializes the examl environment  
   */ 
  void initializeExaml(const CommandLine &cl); 
  void initializeExamlAlt(const CommandLine &cl ); 
  /** 
      @brief a printing function that allows to print properly 
   */ 
  static void blockingPrint( Communicator &comm,std::string ss );  

  std::string printLoadBalance(const TreeAln& traln, nat numRun, nat numCoupled ) const ; 

  /** 
      @brief indicates whether the chain belongs to a process
   */ 
  bool isMyChain(nat runid, nat chainIndex) const ; 
  bool isMyChain(nat gRank, nat runid, nat chainIndex )const;
  /** 
      @brief indicates whether the run belongs to a process
   */   
  bool isMyRun(nat runid) const ; 
  bool isMyRun(nat gRank, nat runid) const; 
  /** 
      @brief a generic exit function  
   */ 
  static void genericExit(int code); 
  friend std::ostream& operator<<(std::ostream& out, const ParallelSetup &pl);

  Communicator& getChainComm() { return _chainComm; }
  Communicator& getGlobalComm() {return _globalComm; }
  Communicator& getRunComm() {return _runComm; }

  IncompleteMesh getMesh()const {return _mesh;}
  
  int getRankInRun( std::array<nat,3> coords) const ; 

  void globalBarrier()  ; 

  bool globalBroadcast( bool val,int root ) const; 

  std::array<nat,3> getMyCoordinates() const; 
  
private: 			// METHODS
  auto serializeAllChains(std::vector<CoupledChains> &runs, CommFlag commFlags) const
    -> std::tuple<nat,std::vector<char>> ; 
  auto serializeSwapMatrices(std::vector<CoupledChains>& runs, CommFlag &commFlags ) const 
    -> std::tuple<nat, std::vector<char>> ; 

private:   			// ATTRIBUTES
  nat _runsParallel;
  nat _chainsParallel; 

  Communicator _globalComm; 
  Communicator _runComm; 
  Communicator _chainComm; 

  std::unordered_map<int,Communicator> _commToChains; // intercomms 

  IncompleteMesh _mesh; 
  
  Communicator _chainLeaderComm; 
}; 


#endif
