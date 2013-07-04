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
  CoupledChains(randCtr_t seed, int runNum, string workingdir, int numCoupled,  vector<Chain>& _chains ); 

  void seedChains(); 

  void initializeChains(vector<shared_ptr<TreeAln> > trees, const vector<unique_ptr<AbstractProposal> > &proposals, 
			const vector<shared_ptr<AbstractParameter> > &vars, shared_ptr<LikelihoodEvaluator> eval); 

  ~CoupledChains(); 

  /** @brief run for a given number of generations */
  void run(int numGen); 

  void printSwapInfo();

  void chainInfo(); 

  /** @brief Execute a portion of one run. */
  void executePart(int gensToRun); 

  
  void setPrintFreq(nat t){printFreq = t; }

  void setSwapInterval(nat i) {swapInterval = i; }
  void setSamplingFreq(nat i) {samplingFreq = i; }
  void setHeatIncrement(double temp ) { heatIncrement = temp ; } 
  void setTuneHeat(bool bla){tuneHeat = bla ; }  
  void setTemperature(double temp ){heatIncrement = temp;  } 
  vector<Chain>& getChains() {return chains; } 
  int getNumberOfChains(){return chains.size();}
  void enableHeatTuning(int freq ) { tuneHeat = true; tuneFreq = freq; }  
  void printNexusTreeFileStart(Chain &chain, FILE *fh  );   
  FILE* getTopoFile(){return topoFile; }
  FILE* getParamFile(){return paramFile; }  
  void setRunName(string a) {runname = a;  }

private: 			// METHODS

/**
   @brief attempt a MC3 switch of chains 
   @param chains -- the pointer to the first chain struct in the chain-array.      
   NOTICE does not work with priors yet 
*/
  void switchChainState(); 

  /** @brief tunes the temperature parameter */  
  void tuneTemperature();


private: 			// ATTRIBUTES

  vector<Chain> chains; 
  vector<SuccessCounter*> swapInfo;  

  double heatIncrement; 
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
  
  nat  numCoupled; 
}; 

#endif

