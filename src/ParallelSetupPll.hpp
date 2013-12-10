#ifndef _PARALLEL_SETUP_PLL_HPP
#define _PARALLEL_SETUP_PLL_HPP

#include <string>
#include <cassert>
#include <vector>
#include <iosfwd>
#include <sstream>

#include "common.h"
#include "CoupledChains.hpp"
#include "CommFlag.hpp"

class TreeAln; 

class ParallelSetup
{
public: 
  static void genericExit(int code)
  {
    exit(code); 
  }

  bool isMyChain(nat a,nat b) const
  {
    return true; 
  }

  bool isMyRun(nat a) const 
  {
    return true; 
  }

  bool isChainLeader() const
  {
    return true; 
  }

  bool isRunLeader() const 
  {
    return true; 
  }

  bool isGlobalMaster() const 
  {
    return true; 
  }


  void globalBarrier() 
  {
  }

  bool globalBroadcast( bool val,int root ) const
  {
    return val; 
  } 
  
  void printLoadBalance(const TreeAln& traln, nat numRun, nat numCoupled ) const 
  {
    
  } 


  friend std::ostream& operator<<(std::ostream& out, const ParallelSetup &pl)
  {
    return out; 
  }

  void synchronizeChainsAtMaster( std::vector<CoupledChains>& runs, CommFlag commFlags) const
  {
    if( ( commFlags & CommFlag::TREE )  != CommFlag::NOTHING)
      {
	// just for comparability to examl 
	for(auto &run : runs)
	  {
	    for(auto &chain : run.getChains())
	      {
		auto &&oss = std::ostringstream{};
		chain.serializeConditionally(oss, commFlags );
		auto &&iss =  std::istringstream{oss.str()}; 
		chain.deserializeConditionally(iss,  commFlags);
	      }
	  }
      }
  }

  std::string sendRecvChain(CoupledChains &c, nat cAIndex, nat cBIndex, std::string aSer, CommFlag &flags)
  {
    assert(0); 
    return ""; 
  }

  
  static void initialize(int argc, char **argv)
  {
    
  }


  static void finalize()
  {
  }


}; 
#endif
