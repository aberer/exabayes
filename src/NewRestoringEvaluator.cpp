#include "NewRestoringEvaluator.hpp"
#include <cassert>


NewRestoringEvaluator::NewRestoringEvaluator()
{
  assert(0); 
} 

double NewRestoringEvaluator::evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain& root,  const std::vector<nat>& partitions, bool fullTraversal)  
{
  assert(0); 
} 


void NewRestoringEvaluator::evalSubtree( TreeAln &traln, const BranchPlain &evalBranch)     
{
  assert(0); 
} 


double NewRestoringEvaluator::evaluate(TreeAln &traln, const BranchPlain &evalBranch,  bool fullTraversal )  
{
  assert(0); 
} 


void NewRestoringEvaluator::imprint(const TreeAln &traln) 
{
  assert(0); 
} 


void NewRestoringEvaluator::resetToImprinted(TreeAln &traln) 
{
  assert(0); 
}    


void NewRestoringEvaluator::resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions) 
{
  assert(0); 
} 


std::unique_ptr<LikelihoodEvaluator> NewRestoringEvaluator::clone() const 
{
  return std::unique_ptr<LikelihoodEvaluator>{new NewRestoringEvaluator(*this)};
}
