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


void RateHetParameter::printSample(std::ostream& fileHandle, const TreeAln &traln) const 
{
  auto content = extractParameter(traln); 
  // TODO numeric limits
  fileHandle << content.values[0] ; 
}

 
void RateHetParameter::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const 
{
  fileHandle << "alpha{"  ;
  bool isFirst = true; 
  for(auto &p : partitions) 
    {
      fileHandle << (isFirst ? "" : "," ) << p ; 
      isFirst = false; 
    }
  fileHandle<< "}" ; 
} 
