/**
   @brief represents a chain 

   notice: the chain currently is in an intermediate stage to full c++
   conversion. That's why we have two public regions. 
 */ 

#ifndef _CHAIN_H
#define  _CHAIN_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include "PriorBelief.hpp"
#include "LikelihoodEvaluator.hpp"
#include "AbstractProposal.hpp"

#include "TopologyFile.hpp"
#include "ParameterFile.hpp"
#include "Checkpointable.hpp"



class TreeAln; 
class AbstractProposal; 


class Chain : public Checkpointable
{
public: 
  ////////////////
  // LIFE CYCLE //
  ////////////////
  Chain(randKey_t seed, std::shared_ptr<TreeAln> _traln, const std::vector<std::unique_ptr<AbstractProposal> > &_proposals, std::unique_ptr<LikelihoodEvaluator> eval);   
  Chain( Chain&& rhs) ; 
  Chain& operator=(Chain rhs) ; 
  
  /** 
      @brief set seet for the chain specific RNG 
   */ 
  void reseed(randKey_t c) { chainRand = Randomness(c); }  
  /** 
      @brief apply saved parameter contents to the tree structure
      @param eval indicates whether an evaluation should be performed after resuming
      @param checkLnl a hack: disable check for the exact same likelihood. Reason for this is resuming a run from a checkpoint with ExaML. It is just extremely hard to get the exact same likelihood 
   */   
  void resume(bool eval, bool checkLnl)  ; 
  /**
     @brief saves the all parameters that are integrated over,
     s.t. the tree can be used by another chain     
     @param paramsOnly indicates whether the likelihood and prior density should be saved as well   
   */ 
  void suspend(bool paramsOnly)  ; 
  /** 
      @brief proceed by one generation 
   */ 
  void step();
  /** 
      @brief gets the proposals of this chain
   */ 
  const std::vector<AbstractProposal*> getProposalView() const  ; 
  /** 
      @brief switch state with chain rhs
   */ 
  void switchState(Chain &rhs);
  /** 
      @brief add a representation of the chain to the stream    
   */ 
  std::ostream& addChainInfo(std::ostream &out) const; 
  /** 
      @brief print success rates of proposals 
   */ 
  void printProposalState(std::ostream& out ) const ; 
  /** 
      @brief extracts the variables of a chain into a sorted array      
   */ 
  const std::vector<AbstractParameter*> extractVariables() const ; 

  // getters and setters 
  double getLnLikelihood() const {return   tralnPtr->getTr()->likelihood;} 
  double getLnPrior() const {return prior.getLnPrior(); }
  double getBestState() const {return bestState; }
  LikelihoodEvaluator* getEvaluator() {return evaluator.get(); }
  const TreeAln& getTraln() const { return *tralnPtr; }
  TreeAln& getTraln()  { return *tralnPtr; }
  Randomness& getChainRand(){return chainRand;}
  std::shared_ptr<TreeAln> getTralnPtr()   {return tralnPtr; } // BAD 
  double getChainHeat(); 
  void setDeltaT(double dt){deltaT = dt; }
  int getCouplingId() const {return couplingId; }
  void setCouplingId(int id) {couplingId = id; }
  void setTuneFreuqency(nat _tuneFreq ) {tuneFrequency = _tuneFreq; }
  void setHeatIncrement(nat cplId) {couplingId = cplId  ;}
  void setRunId(nat id) {runid = id; }
  double getDeltaT() {return deltaT; }
  int getGeneration() const {return currentGeneration; }
  double getLikelihood() const {return likelihood; }
  const PriorBelief& getPrior() const  {return prior; } 
  void sample( const TopologyFile &tFile, const ParameterFile &pFile  ) const ; 

  virtual void readFromCheckpoint( std::ifstream &in ) ; 
  virtual void writeToCheckpoint( std::ofstream &out) ;   

private : 			// METHODS 
  void debug_printAccRejc(AbstractProposal* prob, bool accepted, double lnl, double lnPr, double hastings ) ;
  AbstractProposal* drawProposalFunction();

  void printArrayStart(); 


private: 			// ATTRIBUTES
  std::shared_ptr<TreeAln> tralnPtr;  
  double deltaT; 		// this is the global heat parameter that defines the heat increments  
  int runid; 
  int tuneFrequency; 		// TODO should be have per-proposal tuning?   
  double hastings;  		// logged!
  int currentGeneration;     	
  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  int couplingId;					     // CHECKPOINTED 
  std::vector<std::unique_ptr<AbstractProposal> > proposals; 
  Randomness chainRand;
  double relWeightSum; 	// sum of all relative weights
  PriorBelief prior; 
  double bestState; 
  std::unique_ptr<LikelihoodEvaluator> evaluator;   
  
  // suspending and resuming the chain   
  double likelihood;
  double lnPr; 	

  std::unordered_map<nat, ParameterContent> savedContent; // maps parameter-id to content
  

  // friends 
  friend std::ostream& operator<<(std::ostream& out, const Chain &rhs); 
}; 




#endif
