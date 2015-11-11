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

#include "ProposalSet.hpp"
#include "PriorBelief.hpp"
#include "eval/LikelihoodEvaluator.hpp"
#include "proposals/AbstractProposal.hpp"

#include "file/TopologyFile.hpp"
#include "file/ParameterFile.hpp"
#include "Serializable.hpp"

#include "ParallelSetup.hpp"

class TreeAln; 
class AbstractProposal; 


class Chain : public Serializable
{
public: 
  Chain(randKey_t seed, TreeAln _traln, const std::vector<std::unique_ptr<AbstractProposal> > &_proposals, 
	std::vector<ProposalSet> proposalSets, LikelihoodEvaluator eval, bool isDryRun);   
  Chain( const Chain& rhs)   ; 
  Chain( Chain&& rhs) ; 
  Chain& operator=( Chain rhs) ; 

  /** 
      @brief apply saved parameter contents to the tree structure
      @param eval indicates whether an evaluation should be performed after resuming
      @param checkLnl a hack: disable check for the exact same likelihood. Reason for this is resuming a run from a checkpoint with ExaML. It is just extremely hard to get the exact same likelihood 
   */   
  void resume()  ; 
  /**
     @brief saves the all parameters that are integrated over,
     s.t. the tree can be used by another chain     
     @param paramsOnly indicates whether the likelihood and prior density should be saved as well   
   */ 
  void suspend()  ; 
  /** 
      @brief proceed by one generation 
   */ 
  void step();
  /** 
      @brief gets the proposals of this chain
   */ 
  const std::vector<AbstractProposal*> getProposalView() const  ; 
  /** 
      @brief add a representation of the chain to the stream    
   */ 
  std::ostream& addChainInfo(std::ostream &out) const; 
  /** 
      @brief extracts the variables of a chain into a sorted array      
   */ 
  std::vector<AbstractParameter*> extractParameters() const ; 
  std::vector<AbstractParameter*> extractSortedParameters() const ; 
  /** 
      @brief take a sample from the chain 
   */ 
  void sample( std::unordered_map<nat,TopologyFile> &tFile, ParameterFile &pFile  ) const ; 
  /** 
      @brief deserialize the input string based on the flags 
   */ 
  void deserializeConditionally(std::string str, CommFlag commFlags); 
  /** 
      @brief serializes the chain into a string based on the flags
   */ 
  std::string serializeConditionally( CommFlag commFlags) const ; 

  void reinitPrior(){_prior.initialize(_traln, extractParameters());}

  // getters and setters 
  double getBestState() const {return _bestState; }
  LikelihoodEvaluator& getEvaluator() {return _evaluator; }
  const TreeAln& getTralnHandle() const { return _traln; }
  TreeAln& getTralnHandle()  { return _traln; }
  Randomness& getChainRand(){return _chainRand;}
  double getChainHeat() const; 
  void setDeltaT(double dt){_deltaT = dt; }
  int getCouplingId() const {return _couplingId; }
  void setCouplingId(int id) {_couplingId = id; }
  void setTuneFreuqency(nat _tuneFreq ) {_tuneFrequency = _tuneFreq; }
  void setHeatIncrement(nat cplId) {_couplingId = cplId  ;}
  void setRunId(nat id) {_runid = id; }
  double getDeltaT() {return _deltaT; }
  int getGeneration() const {return _currentGeneration; }
  double getLikelihood() const {return _traln.getTrHandle().likelihood; }
  void setLikelihood(double lnl) { _traln.getTrHandle().likelihood = lnl; } 
  void setLnPr(double lnPr) { _lnPr = lnPr;  }
  double getLnPr() const {return _lnPr; }
  const PriorBelief& getPrior() const  {return _prior; } 
  void updateProposalWeights(); 
 
  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const ;   
  /**
     @brief the chains determines the virtual root needed by the next
     generation
   */
  BranchPlain peekNextVirtualRoot(TreeAln &traln, Randomness rand) ; 

  friend void swap(Chain &lhs, Chain &rhs); 
  friend void swapHeatAndProposals(Chain &chainA, Chain& chainB) ; 

private : 			// METHODS 
  AbstractProposal& drawProposalFunction(Randomness &rand);
  ProposalSet& drawProposalSet( Randomness &rand); 
  void printArrayStart(); 
  void initProposalsFromStream(std::istream& in);

  void stepSingleProposal(); 
  void stepSetProposal();

private: 			// ATTRIBUTES
  TreeAln _traln; 
  double _deltaT; 		// this is the global heat parameter that defines the heat increments  
  int _runid; 
  int _tuneFrequency; 		// TODO should be have per-proposal tuning?   
  double _hastings;  		// logged!
  int _currentGeneration;     	
  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  nat _couplingId;					     // CHECKPOINTED 
  std::vector<std::unique_ptr<AbstractProposal> > _proposals; 
  std::vector<ProposalSet> _proposalSets; 
  Randomness _chainRand;
  double _relWeightSumSingle;
  double _relWeightSumSets;   
  PriorBelief _prior; 
  double _bestState; 
  LikelihoodEvaluator _evaluator;   
  
  // suspending and resuming the chain   
  double _lnPr; 	

  // friends 
  friend std::ostream& operator<<(std::ostream& out, const Chain &rhs); 
}; 

#endif
