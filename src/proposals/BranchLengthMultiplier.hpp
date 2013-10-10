#ifndef _BRANCHLENGTHSMULTIPLIER_H
#define _BRANCHLENGTHSMULTIPLIER_H

#include "AbstractProposal.hpp"

#include <limits>


class BranchLengthMultiplier : public AbstractProposal
{
public: 
  BranchLengthMultiplier(  double multiplier); 
  virtual ~BranchLengthMultiplier(){}

  virtual BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const ; 
  // { return proposeBranch(traln, rand);  } 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln) ; 

  virtual void autotune();

  virtual std::vector<nat> getInvalidatedNodes(const TreeAln &traln ) const  {return {}; }

  virtual AbstractProposal* clone() const;  
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand) ;

  virtual BranchPlain proposeBranch(const TreeAln &traln, Randomness &rand) const ;   

  virtual void readFromCheckpointCore(std::istream &in); 
  virtual void writeToCheckpointCore(std::ostream &out) const; 

protected: 
  double multiplier;  
  BranchLength savedBranch;   

}; 


#endif
