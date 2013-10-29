#ifndef _AAMODEL_JUMP_H
#define _AAMODEL_JUMP_H

#include <vector>

#include "Branch.hpp"
#include "AbstractProposal.hpp"
// #include "axml.h"
#include "ProtModel.hpp"


class AminoModelJump : public AbstractProposal
{
public: 

  AminoModelJump(std::vector<ProtModel> matrices); 

  virtual BranchPlain determinePrimeBranch(const TreeAln& traln, Randomness &rand) const  { return BranchPlain(); } 

  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out) const  { } //disabled

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand , LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln); 
  virtual void autotune()  ;
  virtual AbstractProposal* clone() const ;  
  virtual std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const  {return {}; }

  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0)); }

private: 
  std::vector<ProtModel> matrices; 
}; 

#endif
