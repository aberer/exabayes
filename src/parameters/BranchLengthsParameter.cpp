#include "BranchLengthsParameter.hpp"

void BranchLengthsParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &b : savedContent.branches)
    traln.setBranch(b); 
}


ParameterContent BranchLengthsParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  ParameterContent result; 
  result.branches = traln.extractBranches(); 
  return result; 
}   
