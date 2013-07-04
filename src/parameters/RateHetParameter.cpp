#include "RateHetParameter.hpp"


void RateHetParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  for(auto &m : partitions)
    traln.setAlpha(content.values[0], m); 
}
 

ParameterContent RateHetParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = {traln.getAlpha(partitions[0])} ;
  return result;  
}   
