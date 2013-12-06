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
   them hot
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

  /** 
      @brief indicates whether this run is executed by the process
   */ 
  bool isMyRun (ParallelSetup &pl) const ; 

  /** 
      @brief Execute a portion of one run. 
  */
  void executePart(nat startGen, nat numGen,  ParallelSetup &pl);   
  void setSamplingFreq(nat i) {_samplingFreq = i; }
  void setHeatIncrement(double temp ) { _heatIncrement = temp ; } 
  void setTemperature(double temp ){_heatIncrement = temp;  } 
  vector<Chain>& getChains() {return _chains; } 
  nat getRunid()  const {return _runid; }
  const vector<Chain>& getChains() const {return _chains; }
  int getNumberOfChains(){return _chains.size();}
  void setNumSwapsPerGen(double s){_numSwapsPerGen = s; }
  void setRunName(string a) {_runname = a;  }
  void initializeOutputFiles(bool isDryRun)  ; 
  SwapMatrix getSwapInfo() const {return _swapInfo; }
  void addToSwapMatrix(SwapMatrix toAdd){ _swapInfo = _swapInfo + toAdd;  }
  const Randomness& getRandomness() const {return _rand; }
  std::vector<std::string> getAllFileNames() const ; 
  void finalizeOutputFiles(const ParallelSetup &pl); 

  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const ;   

  void regenerateOutputFiles(std::string _workdir, std::string prevId) ; 
  
private: 			// METHODS
  /**
     @brief attempt to swap two _chains
  */ 
  void attemptSwap( ParallelSetup &pl); 

private: 			// ATTRIBUTES
  std::vector<Chain> _chains; 
  SwapMatrix _swapInfo; 
  double _heatIncrement; 	
  Randomness _rand; 
  int _runid; 
  int _samplingFreq; 
  string _runname; 
  string _workdir; 

  // order is coupled to the heat id  
  std::unordered_map<nat,TopologyFile> _paramId2TopFile; 
  std::vector<ParameterFile> _pFile; 

  double _numSwapsPerGen; 
}; 

#endif

