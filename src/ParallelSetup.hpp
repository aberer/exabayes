/** 
    @file ParallelSetup.hpp
    @brief handles all high-level parallel communication    
 */ 

#ifndef _PARALLEL_SETUP_H
#define _PARALLEL_SETUP_H

#include <unordered_map>

#include "config/CommandLine.hpp"

class SampleMaster;
class CoupledChains; 

#if HAVE_PLL == 0
extern MPI_Comm comm; 		// the examl communicator 
extern int processID; 		// examl rank 
extern int processes; 		// examl comm size 
#endif

#include <iostream>



/** 
    @brief indicates how much communication is necessary     
 */ 
enum class CommFlag : int
{
  PrintStat    =  1 , 		// stats for printing  
    Proposals  =  2 , 		// all proposal data  
    Tree       =  4 ,		// the tree state 
    Swap       =  8    		// swap matrix 
}; 


inline CommFlag operator|(CommFlag a, CommFlag b)
{
  return static_cast<CommFlag>(static_cast<int>(a) | static_cast<int>(b)); 
}


inline bool operator&(CommFlag a, CommFlag b)
{
  return (static_cast<int>(a) & static_cast<int>(b)) != 0 ; 
}


class ParallelSetup
{
public: 
  ParallelSetup(int argc, char **argv);

  /** 
      @brief finalize the parallel environment 
   */ 
  void finalize(); 
  /** 
      @brief gets the number of procesess working on this chain batch 
   */   
  nat getProcessesPerChainBatch() const  { assert(globalSize % (runsParallel * chainsParallel) == 0 ); return globalSize / (runsParallel * chainsParallel); }
  /** 
      @briet gets the id of the batch of chains these processes are working on   
   */ 
  nat getChainBatch() const {return getChainBatch(globalRank);  }
  nat getChainBatch(nat gRank) const {return (gRank / getProcessesPerChainBatch()) % chainsParallel; }  
  /** 
      @brief gets the batch of runs assigned to this proess
   */ 
  nat getRunBatch() const {return getRunBatch(globalRank); } 
  nat getRunBatch(nat gRank) const {return gRank / getProcessesPerRunBatch(); }
  /** 
      @brief gets the number of runs executed in parallel  
   */ 
  nat getRunsParallel() const {return runsParallel ; }
  /** 
      @brief gets the number of (coupled) chains executed in paralel (in addition to run-level parallelism)
   */ 
  nat getChainsParallel() const {return chainsParallel; }
  /** 
      @brief gets number of proceses per batch of (independent) runs
   */ 
  nat getProcessesPerRunBatch() const {assert(globalSize % (runsParallel * chainsParallel) == 0 ); return globalSize / runsParallel; }
  /** 
      @brief gets the rank with respect to the batch of runs 
   */ 
  nat getRankInRunBatch() const {return getRankInRunBatch(globalRank);  }
  nat getRankInRunBatch(nat gRank ) const {return gRank % getProcessesPerRunBatch(); }
  /**
     @brief gets the rank of the process relative to all processes working on this batch of chains 
   */ 
  nat getRankInChainBatch() const {return getRankInChainBatch(globalRank); } 
  nat getRankInChainBatch(nat gRank)const {return gRank % getProcessesPerChainBatch(); }  
  /** 
      @brief indicates whether the process should conduct output operatons (for his run) 
   */ 
  bool isRunLeader() const {  return isRunLeader(globalRank)  ; }  
  bool isRunLeader(nat gRank) const {  return getRankInRunBatch(gRank) == 0 ; }  
  /**
     @brief indicates whether a process is the leader in a batch of chains 
   */ 
  bool isChainLeader() const {return getRankInChainBatch() == 0; }
  /** 
      @brief indicates whether the process is the process to conduct global logging output
   */ 
  bool isGlobalMaster() const {return globalRank == 0 ; }
  /** 
      @brief sends a serialized chain representation to the other involved process  
      @param run relevant run  
      @param otherChainIndex the index of the chain that is involved in the swapping attempt   
      @return serialized representation  of the chain we swap with 
   */ 
  std::string
  sendRecvChain(const CoupledChains& run, nat myIndex, nat otherChainIndex  , std::string myChainSer, CommFlag flags )  ; 

  /** 
      @brief synchronize all core chain information at the master node
      @notice this only concerns the global masterprocess in a
      parallel environtment.  
  */ 
  void synchronizeChainsAtMaster( std::vector<CoupledChains>& runs, CommFlag commFlags) const;  

#if HAVE_PLL == 0
  /** 
      @brief initializes the examl environment  
   */ 
  void initializeExaml(const CommandLine &cl); 
  /** 
      @brief a printing function that allows to print properly 
   */ 
  static void blockingPrint(const MPI::Comm &comm,std::string ss );  
#endif

  void printLoadBalance(const TreeAln& traln ) const ; 

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

private:   
  nat runsParallel;
  nat chainsParallel; 
  nat globalRank;
  nat globalSize;  

#if HAVE_PLL == 0
  MPI::Intracomm chainLeaderComm; /// contains only the root processes of each chain  
  MPI::Intracomm chainComm; 
  MPI::Intracomm runComm; 
  std::unordered_map<int,MPI::Intercomm> commToChains;   
#endif
}; 



#endif
