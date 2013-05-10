/**
   @brief represents a chain 

   notice: the chain currently is in an intermediate stage to full c++
   conversion. That's why we have two public regions. 
 */ 



#ifndef _CHAIN_H
#define  _CHAIN_H
#include <vector>

#include "ConfigReader.hpp"
#include "Category.hpp"
#include "PriorBelief.hpp"




class AbstractProposal; 

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

  /* substitution parameter */
  perPartitionInfo *infoPerPart; 

} paramDump; 


class Chain
{
  // LEGACY stuff 
public: 

  // TODO shared pointer!
  TreeAln *traln; 

  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  int currentGeneration;   

  double penaltyFactor; //change to the probability of picking a proposal  

  /* saves the entire space in the parameter space  */
  paramDump dump;

  double hastings;/// the proposal ratio 

  vector<Category> proposalCategories; // proposals that we implemented using the new framework 

  // CORRECT part 
public: 
  Chain(randKey_t seed, int id, int _runid, TreeAln* _traln, const PriorBelief &_prior, const vector<Category> &propCats) ;
  Chain& operator=(Chain& rhs); 

  // getters and setters 
  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }
  void setRestorer(LnlRestorer *rest){restorer = rest ; }
  LnlRestorer* getRestorer(){return restorer; }
  vector<Category>& getProposalCategories(){return proposalCategories;}
  
  /** @brief draws a proposal function. */ 
  AbstractProposal* drawProposalFunction( ); 

  /** @brief Saves all relevan information from the tree into the chain Chain. */ 
  void saveTreeStateToChain(); 

  /** @brief Applies the Chain of the chain to its tree. */ 
  void applyChainStateToTree(); 
  
  /** @brief Execute one generation of the chain. */
  void step();

  int getGeneration(){return currentGeneration; }

  Randomness* getChainRand(){return chainRand;}

  void printNexusTreeFileStart( FILE *fh  ); 
  void printParams(FILE *fh); 
  void printParamFileStart(FILE *fh); 
  void finalizeOutputFiles(FILE *fh); 
  void printSample(FILE *topofile, FILE *paramFile); 
  void printTopology(FILE *fh); 

  void clarifyOwnership(); 
  
private : 
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  Randomness *chainRand; 
  LnlRestorer *restorer; 
  int runid; 
  PriorBelief prior; 
  
  int tuneFrequency; 		// TODO should be have per-proposal tuning?   

  void debug_printAccRejc(AbstractProposal *prob, bool accepted, double lnl ) ; 
  void initParamDump(); 
}; 


#endif
