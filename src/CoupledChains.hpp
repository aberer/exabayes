/**
   @file CoupledChains.hpp 

   represents a run (consisting of a number of coupled chains) 
   
 */


#ifndef _COUPLED_CHAINS_H
#define _COUPLED_CHAINS_H


#include "TreeAln.hpp"
#include "Randomness.hpp"

class SuccessCtr; 



/**
   @brief represents some coupled chains, one of them cold, many of
   them Hot
 */ 


class CoupledChains
{
public: 
  CoupledChains(int seed, int numCoupled, vector<TreeAln*> trees, int runid, initParamStruct *initParam);   

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

  Chain* getChain(int i) {return chains[i]; }
  int getNumberOfChains(){return chains.size();}
  

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
  vector<SuccessCtr*> swapInfo;  

  double temperature; 
  
  // the mcmc specific randomness used to initialize chains and do the swapping 
  Randomness rand; 

  int runid; 

  
}; 



#endif

