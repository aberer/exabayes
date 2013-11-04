/**
   @file CoupledChains.hpp 

   represents a run (consisting of a number of coupled chains) 
   
 */


#ifndef _COUPLED_CHAINS_H
#define _COUPLED_CHAINS_H

#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "config/BlockRunParameters.hpp"
#include "Chain.hpp"
#include "SuccessCounter.hpp"
#include "file/TopologyFile.hpp"
#include "file/ParameterFile.hpp"
#include "ParallelSetup.hpp"
#include "SwapMatrix.hpp"


/**
   @brief represents some coupled chains, one of them cold, many of
   them Hot
 */ 


class CoupledChains : public Serializable
{
public: 
  ////////////////
  // LIFE CYCLE //
  ////////////////
  CoupledChains(Randomness rand, int runNum, string workingdir, std::string runname, int numCoupled,  vector<Chain>& _chains ); 
  CoupledChains(CoupledChains&& rhs); 
  CoupledChains& operator=(CoupledChains rhs); 

  /**
     @brief run for a given number of generations
  */
  void run(int numGen); 
#if 0 
  void chainInfo(); 
#endif

  /** 
      @brief indicates whether this run is executed by the process
   */ 
  bool isMyRun (ParallelSetup &pl) const ; 

  /** 
      @brief Execute a portion of one run. 
  */
  void executePart(nat startGen, nat numGen,  ParallelSetup &pl);   
  void setPrintFreq(nat t){printFreq = t; }
  void setSwapInterval(nat i) {swapInterval = i; }
  void setSamplingFreq(nat i) {samplingFreq = i; }
  void setHeatIncrement(double temp ) { heatIncrement = temp ; } 
  void setTemperature(double temp ){heatIncrement = temp;  } 
  vector<Chain>& getChains() {return chains; } 
  nat getRunid()  const {return runid; }
  const vector<Chain>& getChains() const {return chains; }
  int getNumberOfChains(){return chains.size();}
  void setNumSwaps(nat ns) {numSwaps = ns; }
  void setRunName(string a) {runname = a;  }
  void initializeOutputFiles(bool isDryRun)  ; 
  SwapMatrix getSwapInfo() const {return swapInfo; }
  void addToSwapMatrix(SwapMatrix toAdd){ swapInfo = swapInfo + toAdd;  }
  const Randomness& getRandomness() const {return rand; }

  std::vector<std::string> getAllFileNames() const ; 

  void finalizeOutputFiles()  ; 

  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const ;   

  void regenerateOutputFiles(std::string workdir, std::string prevId) ; 
  
private: 			// METHODS
  /**
     @brief attempt to swap two chains
  */ 
  void attemptSwap( ParallelSetup &pl); 

private: 			// ATTRIBUTES
  vector<Chain> chains; 
  SwapMatrix swapInfo; 
  double heatIncrement; 	// not checkpointed 
  Randomness rand; 
  int runid; 
  int tuneFreq; 
  int printFreq; 
  int swapInterval; 
  int samplingFreq; 
  string runname; 
  string workdir; 

  // order is coupled to the heat id  
  std::unordered_map<nat,TopologyFile> paramId2TopFile; 
  std::vector<ParameterFile> pFile; 

  nat numSwaps; 
}; 

#endif

