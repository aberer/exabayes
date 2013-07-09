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
#include "LikelihoodEvaluator.hpp"
#include "AbstractProposal.hpp"

#include "TopologyFile.hpp"
#include "ParameterFile.hpp"


class TreeAln; 
class AbstractProposal; 


class Chain
{
public: 
  // life cycle related   
  Chain(randKey_t seed, std::shared_ptr<TreeAln> _traln, const std::vector<std::unique_ptr<AbstractProposal> > &_proposals, std::shared_ptr<LikelihoodEvaluator> eval); 
  
  Chain( Chain&& rhs) ; 
  Chain& operator=(Chain rhs) ; 

  void reseed(randKey_t c) { chainRand = Randomness(c); }

  
  void resume(); 

  /** @brief saves the current model state in the parameter objecets that integrate over this */ 
  void suspend(); 
  
  /** @brief draws a proposal function. */ 
  AbstractProposal* drawProposalFunction();

  /** @brief Execute one generation of the chain. */
  void step();

  const std::vector<AbstractProposal*> getProposalView() const  ; 
  void switchState(Chain &rhs);
  std::ostream& addChainInfo(std::ostream &out) const; 
  void printProposalSate(std::ostream& out ) const ; 
  void printProposalState(std::ostream& out ) const ; 
  /** 
      @brief extracts the variables of a chain into a sorted array      
   */ 
  const std::vector<AbstractParameter*> extractVariables() const ; 

  // getters and setters 
  double getLnLikelihood() const {return   traln->getTr()->likelihood;} 
  double getLnPrior() const {return prior.getLnPrior(); }
  double getBestState() const {return bestState; }
  LikelihoodEvaluator& getEvaluator() {return *evaluator; }
  std::shared_ptr<LikelihoodEvaluator> getEvaluatorPtr() {return evaluator; }
  const TreeAln& getTraln() const { return *traln; }
  TreeAln& getTraln()  { return *traln; }
  Randomness& getChainRand(){return chainRand;}
  std::shared_ptr<TreeAln> getTralnPtr()   {return traln; }
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


  void sample( const TopologyFile &tFile, const ParameterFile &pFile  ) const ; 
  // {
  //   tFile.sample( *traln, getGeneration() ); 
  //   pFile.sample( *traln, extractVariables(), getGeneration(), prior.getLnPrior()); 
  // }



private : 			// METHODS 
  void initParamDump(); 	// BAD
  void debug_printAccRejc(AbstractProposal* prob, bool accepted, double lnl, double lnPr ) ;

private: 			// ATTRIBUTES
  std::shared_ptr<TreeAln> traln;  
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  int runid; 
  int tuneFrequency; 		// TODO should be have per-proposal tuning?   
  double hastings;/// the log hastings ratio  
  int currentGeneration;     
  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  std::vector<std::unique_ptr<AbstractProposal> > proposals; 
  Randomness chainRand;   
  double relWeightSum ; 	// sum of all relative weights
  PriorBelief prior; 
  double bestState; 
  std::shared_ptr<LikelihoodEvaluator> evaluator;   
  
  double likelihood; 
  double lnPr; 

  // friends 
  friend std::ostream& operator<<(std::ostream& out, const Chain &rhs); 
}; 




#endif
