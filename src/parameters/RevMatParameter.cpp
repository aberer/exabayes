#include "RevMatParameter.hpp"
#include "axml.h"

void RevMatParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &m : model)
    traln.setRevMat(content.values, m);
} 


ParameterContent RevMatParameter::extractParameter(const TreeAln &traln, std::vector<nat> model)  const
{
  ParameterContent result; 
  result.values = traln.getRevMat(model[0]); 
  return result; 
}   

