#include "BranchLengthsParameter.hpp"

void BranchLengthsParameter::applyParameter(TreeAln& traln, vector<nat> model, ParameterContent &content) const
{
  for(auto &b : savedContent.branches)
    b.applyToTree(traln);
}


ParameterContent& BranchLengthsParameter::extractParameter(const TreeAln &traln, vector<nat> model)  const
{
  assert(0); 
  // TODO assert 

  savedContent.branches = traln.extractBranches();
}   
