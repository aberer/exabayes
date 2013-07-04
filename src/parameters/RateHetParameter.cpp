#include "RateHetParameter.hpp"


void RateHetParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &m : model)
    traln.setAlpha(content.values[0], m); 
}
 

ParameterContent RateHetParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  ParameterContent result; 
  result.values = {traln.getAlpha(model[0])} ;
  return result;  
}   
