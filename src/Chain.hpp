#ifndef _CHAIN_H
#define  _CHAIN_H
#include <vector>

#include "rng.h"

#include "nclConfigReader.h"

using namespace std; 

class TreeAln; 
class LnlRestorer; 
class Topology; 
class Randomness; 

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
  Topology *topology; 

  /* substitution parameter */
  perPartitionInfo *infoPerPart; 

} paramDump; 







class Chain
{
  // LEGACY stuff 

public: 

  // TODO shared pointer!
  TreeAln *traln; 

  // int id;   
  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  int currentGeneration;   
  
  double priorProb; /// the prior probability of the current state   => store in log form  

  double hastings;/// the proposal ratio 

  double penaltyFactor; //change to the probability of picking a proposal  

  proposalFunction **proposals; 
  int numProposals; 
  double *categoryWeights; 

  /* saves the entire space in the parameter space  */
  paramDump dump;

  /* new stuff that we need when having multiple chains  */
  FILE *topologyFile; 
  FILE *outputParamFile; 



  // CORRECT part 
public: 
  Chain(randKey_t seed, int id, int runid, TreeAln* trealns, initParamStruct *initParams)  ; 
  Chain& operator=(Chain& rhs); 

  void setupProposals(initParamStruct *initParam);

  // getters and setters 
  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }
  void setRestorer(LnlRestorer *rest){restorer = rest ; }
  LnlRestorer* getRestorer(){return restorer; }
  
  /** @brief draws a proposal function. */ 
  void drawProposalFunction(proposalFunction **result ); 

  /** @brief Saves all relevan information from the tree into the chain Chain. */ 
  void saveTreeStateToChain(); 


  /** @brief Applies the Chain of the chain to its tree. */ 
  void applyChainStateToTree(); 
  
  /** @brief Execute one generation of the chain. */
  void step();

  Randomness* getChainRand(){return chainRand;}

  void printInfo(const char *format, ...); 
  
private : 
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  Randomness *chainRand; 
  LnlRestorer *restorer; 
  int runid; 
  
  void normalizePropSubCats(); 
  void normalizeCategories(); 
  void printAllProposalWeights(); 
  void debug_printAccRejc(proposalFunction *pf, bool accepted) ; 

  void initParamDump(); 



}; 



#endif
