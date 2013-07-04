#include "BranchLengthsParameter.hpp"

void BranchLengthsParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &b : savedContent.branches)
    b.applyToTree(traln);
}


ParameterContent BranchLengthsParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  assert(0); 
  // TODO assert 

  // savedContent.branches.clear(); 
  // savedContent.branches (traln.extractBranches()) ;
  ParameterContent result; 
  result.branches = traln.extractBranches(); 
  return result; 
}   
