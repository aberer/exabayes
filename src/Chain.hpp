#ifndef _CHAIN_H
#define  _CHAIN_H

#include "rng.h"
#include <vector>

#include "nclConfigReader.h"


using namespace std; 

class TreeAln; 
class LnlRestorer; 
class Topology; 

typedef struct _pfun proposalFunction; 




typedef struct
{
  int modelNum; 
  double alpha ;
 
  /* TODO only works with DNA */
  int numRates; 
  double substRates[6]; 

  int numFreqs; 
  double frequencies[4];

} perPartitionInfo; 		/* relevant info from pInfo  */


typedef struct 
{
  /* topology */
  // topol *topo; 
  Topology *topology; 
  
  /* branch lengths  */
  double *branchLengths; 

  /* substitution parameter */
  perPartitionInfo *infoPerPart; 

} paramDump; 





class Chain
{
  // LEGACY stuff 

public: 

  // TODO shared pointer!
  TreeAln *traln; 

  bool wasAccepted; 	/// for debug only  

  int id;   
  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  int currentGeneration;   
  
  double priorProb; /// the prior probability of the current state   => store in log form  

  LnlRestorer *restorer; 

  /* lnlContainer lnl;  /// the current likelihood of the state */

  double hastings;/// the proposal ratio 

  double penaltyFactor; //change to the probability of picking a proposal  

  proposalFunction **proposals; 
  int numProposals; 
  double *categoryWeights; 

  /* new stuff that we need when having multiple chains  */
  FILE *topologyFile; 
  FILE *outputParamFile; 

  /* RNG */
  randKey_t rKey;
  randCtr_t rCtr;


  /* saves the entire space in the parameter space  */
  paramDump dump;



  // CORRECT part 
public: 
  Chain(randKey_t seed, int id, int runid, TreeAln* trealns, initParamStruct *initParams)  ; 

  
  void setRestorer(LnlRestorer *rest){restorer = rest ; }
  void setupProposals(initParamStruct *initParam);

  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }
  
  /** @brief Execute one generation of a given chain. */
  void step();
  

private : 
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  
  void normalizePropSubCats(); 
  void normalizeCategories(); 
  void printAllProposalWeights(); 

}; 



#endif
