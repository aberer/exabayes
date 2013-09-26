#ifndef _NEW_RESTORING_EVALUATORm
#define _NEW_RESTORING_EVALUATORm

#include "LikelihoodEvaluator.hpp"




class NewRestoringEvaluator : public LikelihoodEvaluator
{
public: 
  NewRestoringEvaluator(); 
  virtual ~NewRestoringEvaluator(){}

  virtual double evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain& root,  const std::vector<nat>& partitions, bool fullTraversal) ; 
  virtual void evalSubtree( TreeAln &traln, const BranchPlain &evalBranch)    ; 
  virtual double evaluate(TreeAln &traln, const BranchPlain &evalBranch,  bool fullTraversal ) ; 
  virtual void imprint(const TreeAln &traln); 
  virtual void resetToImprinted(TreeAln &traln);    
  virtual void resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions); 
  virtual std::unique_ptr<LikelihoodEvaluator> clone() const; 

private: 
  

}; 

#endif
