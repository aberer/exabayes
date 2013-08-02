#ifndef RESTORING_LIKELIHOOD_EVALUATOR
#define RESTORING_LIKELIHOOD_EVALUATOR

#include "LikelihoodEvaluator.hpp"

class RestoringLnlEvaluator : public LikelihoodEvaluator
{
public: 
  RestoringLnlEvaluator(std::shared_ptr<ArrayRestorer> restorer); 
  
  virtual ~RestoringLnlEvaluator(){}; 

  double evaluateNoBack(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 

  virtual double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal ) ;  
  virtual void evalSubtree( TreeAln &traln, const Branch &evalBranch); 
  virtual double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions, bool fullTraversal) ; 


  virtual void resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions) ; 

  virtual void resetToImprinted(TreeAln &traln); 
  virtual void imprint(const TreeAln &traln); 

  virtual std::unique_ptr<LikelihoodEvaluator> clone() const {return std::unique_ptr<LikelihoodEvaluator>(new  RestoringLnlEvaluator(*this));  } 

private: 
  std::shared_ptr<ArrayRestorer> restorer;    
}; 


#endif
