/**
   @file CoupledChains.hpp 

   represents a run (consisting of a number of coupled chains) 
   
 */


#ifndef _COUPLED_CHAINS_H
#define _COUPLED_CHAINS_H

#include <queue>

#include "model/TreeAln.hpp"
#include "math/Randomness.hpp"
#include "config/BlockRunParameters.hpp"
#include "Chain.hpp"
#include "SuccessCounter.hpp"
#include "file/TopologyFile.hpp"
#include "file/ParameterFile.hpp"
#include "SwapMatrix.hpp"
#include "comm/SwapElem.hpp"

class PendingSwap; 
class ParallelSetup;

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
  CoupledChains(CoupledChains&& rhs) = default; 
  CoupledChains(const CoupledChains& rhs) = delete; 
  CoupledChains& operator=(CoupledChains rhs); 

  friend void swap(CoupledChains &lhs, CoupledChains& rhs ); 

  // void deleteMyFiles() const ; 

  /**
     @brief run for a given number of generations
  */
  void run(int numGen); 

  /** 
      @brief Execute a portion of one run. 
  */
  void executePart(nat startGen, nat numGen,  ParallelSetup &pl);   
  void  executePartNew(nat startGen, nat numGen, ParallelSetup& pl); 
  void doStep(nat id, ParallelSetup &pl); 
  void setSamplingFreq(nat i) {_samplingFreq = i; }
  void setHeatIncrement(double temp ) { _heatIncrement = temp ; } 
  void setTemperature(double temp ){_heatIncrement = temp;  } 
  std::vector<Chain>& getChains() {return _chains; } 
  nat getRunid()  const {return _runid; }
  const vector<Chain>& getChains() const {return _chains; }
  int getNumberOfChains(){return _chains.size();}
  void setNumSwapsPerGen(double s){_numSwapsPerGen = s; }
  void setRunName(string a) {_runname = a;  }
  void initializeOutputFiles(bool isDryRun)  ; 
  SwapMatrix getSwapInfo() const {return _swapInfo; }
  void setSwapInfo(SwapMatrix swap){ _swapInfo = swap; }
  void addToSwapMatrix(const SwapMatrix &toAdd){ _swapInfo = _swapInfo + toAdd;  }
  const Randomness& getRandomness() const {return _rand; }
  std::vector<std::string> getAllFileNames() const ; 

  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const ;   

  void regenerateOutputFiles(std::string _workdir, std::string prevId) ; 
  std::list<SwapElem> generateSwapsForBatch(nat startGen, nat numGen) ; 
  
private: 			// METHODS
  /**
     @brief attempt to swap two _chains
  */ 
  bool allMyChainsAreBlocked( const std::vector<bool> &isBlocked, const ParallelSetup& pl ) const ; 
  bool doLocalSwap(ParallelSetup &pl, const SwapElem &theSwap ); 
  bool doSwap( ParallelSetup &pl, const SwapElem& elem ); 

  PendingSwap prepareSwap(ParallelSetup &pl, const SwapElem& theSwap); 

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


