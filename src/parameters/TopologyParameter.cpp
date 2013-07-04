#include "TopologyParameter.hpp"


void TopologyParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  // this is not what we want! topology should be restored! 
  assert(0); 
  for(auto & b : savedContent.branches)
    traln.setBranch(b); 
}


ParameterContent TopologyParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  assert(model.size() == 1 ); 

  ParameterContent result; 
  result.branches = traln.extractBranches(); 
  return result; 
}   
