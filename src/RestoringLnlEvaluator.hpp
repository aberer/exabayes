#ifndef RESTORING_LIKELIHOOD_EVALUATOR
#define RESTORING_LIKELIHOOD_EVALUATOR

#include "LikelihoodEvaluator.hpp"

class RestoringLnlEvaluator : public LikelihoodEvaluator
{
public: 
  RestoringLnlEvaluator(std::shared_ptr<LnlRestorer> restorer); 

  double evaluateNoBack(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 

  virtual double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal ) ;  
  virtual void evalSubtree( TreeAln &traln, const Branch &evalBranch); 
  virtual double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions) ; 

  virtual void resetToImprinted(TreeAln &traln); 
  virtual void imprint(const TreeAln &traln); 

private: 
  std::shared_ptr<LnlRestorer> restorer;    

}; 


#endif
