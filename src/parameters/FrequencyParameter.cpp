#include "FrequencyParameter.hpp"


void FrequencyParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  // std::cout << "proposal applies: " ; for_each(content.values.begin(), content.values.end(), [](double d ) {std::cout << d << "d"; }) ; std::cout << std::endl; 
  
  for(auto &m : partitions)    
    traln.setFrequencies(content.values, m); 
}


ParameterContent FrequencyParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = traln.getFrequencies(partitions[0]); 
  return result; 
}   
