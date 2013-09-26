#include "TopologyParameter.hpp"


void TopologyParameter::applyParameter(TreeAln& traln , const ParameterContent &content) const
{
  traln.unlinkTree();
  for(auto &b : content.topology)
    traln.clipNode(traln.getUnhookedNode(b.getPrimNode()), traln.getUnhookedNode(b.getSecNode())); 
}


ParameterContent TopologyParameter::extractParameter(const TreeAln &traln )  const
{
  auto result = ParameterContent{}; 
  auto params = std::vector<AbstractParameter*>{};
  result.topology = traln.extractBranches(  ); 
  return result; 
}   

