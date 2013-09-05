#include "BranchLengthsParameter.hpp"

void BranchLengthsParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  // auto view = getPrimaryParameterView();
  // assert(view.size() == 1); 

  // TODO so it has gotten to this ... => const correctness! 
  
  for(auto &b : content.branchLengths)
    traln.setBranch(b, const_cast<AbstractParameter*>(dynamic_cast<const AbstractParameter* const >(this))); 
}

ParameterContent BranchLengthsParameter::extractParameter(const TreeAln &traln)  const
{
  auto result = ParameterContent{}; 

  // TODO bad 
  AbstractParameter* param = const_cast<AbstractParameter*>(dynamic_cast<const AbstractParameter* const>(this)); 

  result.branchLengths = traln.extractBranches(param); 

  return result; 
}   

