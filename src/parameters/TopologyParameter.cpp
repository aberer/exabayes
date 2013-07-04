#include "TopologyParameter.hpp"


void TopologyParameter::applyParameter(TreeAln& traln , const ParameterContent &content) const
{
  // this is not what we want! topology should be restored! 
  assert(0); 
  for(auto & b : savedContent.branches)
    traln.setBranch(b); 
}


ParameterContent TopologyParameter::extractParameter(const TreeAln &traln )  const
{
  assert(partitions.size() == 1 ); 

  ParameterContent result; 
  result.branches = traln.extractBranches(); 
  return result; 
}   
