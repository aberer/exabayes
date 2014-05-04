/** 
    @brief our flavour of the statistical nearest neighbor interchange

    in constrast to mb, this proposal will always induce a topological
    change
    
    @notice we could multiply some additional branches
 */ 


#ifndef _STAT_NNI_H
#define _STAT_NNI_H

class Chain; 
#include "data-struct/Path.hpp"
#include "AbstractProposal.hpp"
#include "SprMove.hpp"

class StatNNI : public AbstractProposal
{
public: 
  StatNNI( double multiplier);

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln) ; 

  BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const; 
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0) ); }
  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  
  
  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out) const { } //disabled


  std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const;  

private: 			// METHODS
  void treatOneBranch(nodeptr p, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand); 

private:			// ATTRIBUTES
  double multiplier; 
  SprMove move; 
  bool branchesSaved; 
  std::vector<BranchLengths> savedBls;
};

#endif
