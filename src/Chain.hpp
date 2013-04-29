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
  Chain& operator=(Chain& rhs); 

  
  void setRestorer(LnlRestorer *rest){restorer = rest ; }
  void setupProposals(initParamStruct *initParam);

  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }


  /**   @brief draws a proposal function.

	Notice: this could be extended later, if we decide to make this
	dependent on the previous state.
   
	Furthermore, we must be sure now that category weights and relative
	proposal weights sum up to 1 each. 
   
  */ 
  void drawProposalFunction(proposalFunction **result ); 

  /** @brief Saves all relevan information from the tree into the chain Chain. */ 
  void saveTreeStateToChain(); 


  /** @brief Applies the Chain of the chain to its tree. 

      Notice: you must not simply change the tree pointer. More
      modifications are necessary to do so.

      @param boolean checkLnl -- should we check, if the lnl is the same as
      before? If we applied it the first time, there is no before.
  */ 
  void applyChainStateToTree(); 
  
  /** @brief Execute one generation of a given chain. */
  void step();
  
private : 
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  
  void normalizePropSubCats(); 
  void normalizeCategories(); 
  void printAllProposalWeights(); 

  void initParamDump(); 

}; 



#endif
