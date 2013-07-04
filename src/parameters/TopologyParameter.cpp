#include "TopologyParameter.hpp"


void TopologyParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &b : savedContent.branches)
    b.applyToTree(traln);
}


ParameterContent TopologyParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  assert(model.size() == 1 ); 

  ParameterContent result; 
  result.branches = traln.extractBranches(); 
  return result; 
}   
