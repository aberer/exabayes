/**
   @brief represents a chain 

   notice: the chain currently is in an intermediate stage to full c++
   conversion. That's why we have two public regions. 
 */ 

#ifndef _CHAIN_H
#define  _CHAIN_H

#include <vector>
#include <unordered_set>
#include <memory>

#include "PriorBelief.hpp"
#include "State.hpp"
#include "LikelihoodEvaluator.hpp"
#include "AbstractProposal.hpp"

class TreeAln; 
class AbstractProposal; 

using namespace std; 

class Chain
{
public: 
  // life cycle related   
  Chain(randKey_t seed, shared_ptr<TreeAln> _traln, const vector<unique_ptr<AbstractProposal> > &_proposals, shared_ptr<LikelihoodEvaluator> eval); 
  Chain(  const Chain& rhs)    ; 
  Chain& operator=(Chain &rhs) = delete; 

  void reseed(randKey_t c) { chainRand = Randomness(c); }
  void resume(); 
  void suspend(); 
  
  /** @brief draws a proposal function. */ 
  AbstractProposal* drawProposalFunction();

  /** @brief Execute one generation of the chain. */
  void step();

  const vector<AbstractProposal*> getProposalView() const  ; 
  void printNexusTreeFileStart( FILE *fh  ); 
  void printParams(FILE *fh); 
  void printParamFileStart(FILE *fh); 
  void finalizeOutputFiles(FILE *fh); 
  void printSample(FILE *topofile, FILE *paramFile); 
  void printTopology(FILE *fh); 
  void switchState(Chain &rhs);
  ostream& addChainInfo(ostream &out) const; 
  void printProposalSate(ostream& out ) const ; 
  void printProposalState(ostream& out ) const ; 
  vector<RandomVariable*> extractVariables() const ; 

  // getters and setters 
  double getLnLikelihood() const {return   traln->getTr()->likelihood;} 
  double getLnPrior() const {return prior.getLnPrior(); }
  double getBestState() const {return bestState; }
  LikelihoodEvaluator& getEvaluator() {return *evaluator; }
  shared_ptr<LikelihoodEvaluator> getEvaluatorPtr() {return evaluator; }
  const TreeAln& getTraln() const { return *traln; }
  TreeAln& getTraln()  { return *traln; }
  Randomness& getChainRand(){return chainRand;}
  shared_ptr<TreeAln> getTralnPtr()   {return traln; }
  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }
  int getCouplingId(){return couplingId; }
  void setCouplingId(int id) {couplingId = id; }
  void setTuneFreuqency(nat _tuneFreq ) {tuneFrequency = _tuneFreq; }
  void setHeatIncrement(nat cplId) {couplingId = cplId  ;}
  void setRunId(nat id) {runid = id; }
  double getDeltaT() {return deltaT; }
  int getGeneration() const {return currentGeneration; }
  const PriorBelief& getPrior() const  {return prior; } 

private : 			// METHODS 
  void initParamDump(); 	// BAD
  void debug_printAccRejc(AbstractProposal* prob, bool accepted, double lnl, double lnPr ) ;

private: 			// ATTRIBUTES
  shared_ptr<TreeAln> traln;  
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  int runid; 
  int tuneFrequency; 		// TODO should be have per-proposal tuning?   
  double hastings;/// the log hastings ratio  
  int currentGeneration;     
  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  vector<unique_ptr<AbstractProposal> > proposals; 
  State state; 
  Randomness chainRand;   
  double relWeightSum ; 	// sum of all relative weights
  PriorBelief prior; 
  double bestState; 
  shared_ptr<LikelihoodEvaluator> evaluator;   
  double likelihood; 

  // friends 
  friend ostream& operator<<(ostream& out, const Chain &rhs); 
}; 




#endif
