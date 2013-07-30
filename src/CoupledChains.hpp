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


class CoupledChains : public Checkpointable
{
public: 
  ////////////////
  // LIFE CYCLE //
  ////////////////
  CoupledChains(randCtr_t seed, int runNum, string workingdir, std::string runname, int numCoupled,  vector<Chain>& _chains ); 
  CoupledChains(CoupledChains&& rhs); 
  CoupledChains& operator=(CoupledChains rhs); 

  /**
     @brief run for a given number of generations
  */
  void run(int numGen); 
  void chainInfo(); 

  /** 
      @brief indicates whether this run is executed by the process
   */ 
  bool isMyRun (ParallelSetup &pl) const ; 

  /** 
      @brief Execute a portion of one run. 
  */
  void executePart(nat startGen, nat numGen,  ParallelSetup &pl);   
  void setPrintFreq(nat t){printFreq = t; }
  void seedChains(); 
  void setSwapInterval(nat i) {swapInterval = i; }
  void setSamplingFreq(nat i) {samplingFreq = i; }
  void setHeatIncrement(double temp ) { heatIncrement = temp ; } 
  void setTuneHeat(bool bla){tuneHeat = bla ; }  
  void setTemperature(double temp ){heatIncrement = temp;  } 
  vector<Chain>& getChains() {return chains; } 
  nat getRunid()  const {return runid; }
  const vector<Chain>& getChains() const {return chains; }
  int getNumberOfChains(){return chains.size();}
  void enableHeatTuning(int freq ) { tuneHeat = true; tuneFreq = freq; }  
  void printNexusTreeFileStart(Chain &chain, FILE *fh  );   
  void setRunName(string a) {runname = a;  }
  void initializeOutputFiles()  ; 
  SwapMatrix getSwapInfo() const {return swapInfo; }
  void addToSwapMatrix(SwapMatrix toAdd){ swapInfo = swapInfo + toAdd;  }
  const Randomness& getRandomness() const {return rand; }
  void setRandomness(Randomness _rand){rand = _rand; } 

  void finalizeOutputFiles()  ; 

  virtual void readFromCheckpoint( std::istream &in ) ; 
  virtual void writeToCheckpoint( std::ostream &out) const ;   

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
  bool tuneHeat; 
  int printFreq; 
  int swapInterval; 
  int samplingFreq; 
  string runname; 
  string workdir; 

  // order is coupled to the heat id  
  std::vector<TopologyFile> tFile; 
  std::vector<ParameterFile> pFile; 
}; 

#endif

