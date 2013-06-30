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


/**
   @brief represents some coupled chains, one of them cold, many of
   them Hot
 */ 


class CoupledChains
{
public: 
  CoupledChains(int seed, int runNum, const BlockRunParameters &params, vector<TreeAlnPtr> trees, string workingdir, const vector<ProposalPtr> &proposals, const vector<RandomVariablePtr > &vars, LikelihoodEvaluatorPtr eval);



  /** @brief initializes all trees with a given starting tree */
  void initStartingTree(FILE *fh);
  /** @brief initializes all trees with ONE random tree */ 
  void initRandomOne(int seed); 
  
  ~CoupledChains(); 

  /** @brief run for a given number of generations */
  void run(int numGen); 

  void printSwapInfo();

  void chainInfo(); 

  /** @brief Execute a portion of one run. */
  void executePart(int gensToRun); 

  vector<Chain*> getChains() {return chains; } 
  
  // Chain* getChain(int i) {return chains[i]; }
  int getNumberOfChains(){return chains.size();}

  void enableHeatTuning(int freq )
  {
    tuneHeat = true; 
    tuneFreq = freq; 
  }
  
  void printNexusTreeFileStart(Chain &chain, FILE *fh  ); 
  
  FILE* getTopoFile(){return topoFile; }
  FILE* getParamFile(){return paramFile; }
  
  

private: 

/**
   @brief attempt a MC3 switch of chains 
   @param chains -- the pointer to the first chain struct in the chain-array.      
   NOTICE does not work with priors yet 
*/
  void switchChainState(); 


  /** @brief tunes the temperature parameter */  
  void tuneTemperature();


  vector<Chain*> chains; 
  vector<SuccessCounter*> swapInfo;  

  double temperature; 
  // the mcmc specific randomness used to initialize chains and do the swapping 
  Randomness rand; 

  int runid; 

  int tuneFreq; 
  bool tuneHeat; 
  int printFreq; 
  int swapInterval; 
  int samplingFreq; 
  string runname; 

  // files for sampling the cold chain  
  FILE *topoFile; 
  FILE *paramFile; 
}; 

typedef unique_ptr<CoupledChains> CoupledChainsPtr; 

#endif

