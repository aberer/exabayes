#include "BranchLengthsParameter.hpp"



void BranchLengthsParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  for(auto &b : savedContent.branches)
    traln.setBranch(b); 
}

ParameterContent BranchLengthsParameter::extractParameter(const TreeAln &traln)  const
{
  ParameterContent result; 
  result.branches = traln.extractBranches(); 
  return result; 
}   
