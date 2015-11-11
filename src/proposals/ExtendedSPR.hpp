#ifndef _EXTENDED_SPR_H
#define _EXTENDED_SPR_H

#include "AbstractProposal.hpp"
#include "math/Randomness.hpp"
#include "data-struct/Path.hpp"
#include "SprMove.hpp"


class ExtendedSPR : public AbstractProposal
{
public: 
  ExtendedSPR(  double stopProb, double multiplier); 

  virtual BranchPlain determinePrimeBranch( const TreeAln &traln, Randomness &rand) const ; 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln) ; 
  virtual void autotune() {}	// disabled 
  virtual AbstractProposal* clone() const;  
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return 
      std::make_pair(BranchPlain(0,0),BranchPlain(0,0) );}

  virtual void readFromCheckpointCore(std::istream &in) { } 
  virtual void writeToCheckpointCore(std::ostream &out) const { }  

  virtual std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const; 

protected: 			// METHODS
  /**
     @brief draws a random set of branches in the tree that constitute a
     path
   
     This function employs the eSPR strategy 
   
     @param s -- the result : the first two branches in the stack define
     the root of the pruned subtree
  */
  void drawPathForESPR( TreeAln& traln, Randomness &rand, double stopProp ); 

protected: 			// ATTRIBUTES
  double stopProb; 
  double multiplier; 
  SprMove move; 

  bool branchesSaved; 
  std::vector<BranchLengths> savedBls; 
}; 


#endif
