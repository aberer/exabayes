/**
   @brief represents a chain 

   notice: the chain currently is in an intermediate stage to full c++
   conversion. That's why we have two public regions. 
 */ 



#ifndef _CHAIN_H
#define  _CHAIN_H
#include <vector>

using namespace std; 


#include "PriorBelief.hpp"
#include "State.hpp"

class TreeAln; 
class Topology; 
class AbstractProposal; 


class Chain
{
  // LEGACY stuff 
public: 

  // BAD BAD BAD 
  // TODO shared pointer!
  TreeAln *traln; 



public: 
  Chain(randKey_t seed, int id, int _runid, TreeAln* _traln, const vector< unique_ptr<AbstractProposal> > &_proposals, int _tuneFreq ,const vector<RandomVariable> &variables) ; 
  
  // getters and setters 
  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }
  // vector<Category>& getProposalCategories(){return proposalCategories;}
  int getCouplingId(){return couplingId; }
  void setCouplingId(int id) {couplingId = id; }
  
  /** @brief draws a proposal function. */ 
  unique_ptr<AbstractProposal>& drawProposalFunction();

  /** @brief Saves all relevan information from the tree into the chain Chain. */ 
  void saveTreeStateToChain(); 

  /** @brief Applies the Chain of the chain to its tree. */ 
  void applyChainStateToTree(); 
  
  /** @brief Execute one generation of the chain. */
  void step();

  int getGeneration(){return currentGeneration; }

  const PriorBelief& getPrior() const  {return prior; } 

  Randomness& getChainRand(){return chainRand;}

  void printNexusTreeFileStart( FILE *fh  ); 
  void printParams(FILE *fh); 
  void printParamFileStart(FILE *fh); 
  void finalizeOutputFiles(FILE *fh); 
  void printSample(FILE *topofile, FILE *paramFile); 
  void printTopology(FILE *fh); 

  void switchState(Chain &rhs);
  ostream& addChainInfo(ostream &out); 
  
private : 
  Chain(const Chain& rhs)  ; 
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  int runid; 

  int tuneFrequency; 		// TODO should be have per-proposal tuning?   

  double hastings;/// the log hastings ratio  
  int currentGeneration;   
  
  void debug_printAccRejc(unique_ptr<AbstractProposal> &prob, bool accepted, double lnl, double lnPr ) ;
  void initParamDump(); 

  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!

  vector< unique_ptr<AbstractProposal> > proposals; 

  State state; 

  Randomness chainRand; 

  
  double relWeightSum ; 	// sum of all relative weights

  PriorBelief prior; 

}; 


#endif
