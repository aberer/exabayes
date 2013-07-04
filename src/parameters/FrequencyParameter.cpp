#include "FrequencyParameter.hpp"


void FrequencyParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &m : model)    
    traln.setFrequencies(content.values, m); 
}


ParameterContent FrequencyParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  ParameterContent result; 
  result.values = traln.getFrequencies(model[0]); 
  return result; 
}   
