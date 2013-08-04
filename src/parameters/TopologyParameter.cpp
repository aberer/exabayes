#include "TopologyParameter.hpp"


void TopologyParameter::applyParameter(TreeAln& traln , const ParameterContent &content) const
{
  traln.unlinkTree();
  for(auto &b : content.branches)
    traln.clipNode(traln.getUnhookedNode(b.getPrimNode()), traln.getUnhookedNode(b.getSecNode())// , b.getLength()
		   ); 
}


ParameterContent TopologyParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  std::vector<AbstractParameter*> params;
  result.branches = traln.extractBranches(params); 
  return result; 
}   
