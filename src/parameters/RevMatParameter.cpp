#include <algorithm>
#include <functional>

#include "RevMatParameter.hpp"
#include "axml.h"

void RevMatParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  std::vector<double> tmp = content.values; 
  for_each(tmp.begin(), tmp.end(), [&](double &d){return d /= *(tmp.rbegin()) ;  } );

  for(auto &m : partitions)
    traln.setRevMat(tmp, m);
} 


ParameterContent RevMatParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = traln.getRevMat(partitions[0]); 
  double sum = 0; 
  for(auto &v : result.values )
    sum += v;  
  for_each(result.values.begin(), result.values.end(), [=] (double &d){ d/= sum; }); 
  return result; 
}   



void RevMatParameter::printSample(std::ostream& fileHandle, const TreeAln &traln) const 
{
  auto content = extractParameter(traln);

  bool isFirst = true; 
  for(auto &v : content.values)
    {
      fileHandle << (isFirst ?  "" : "\t" ) << v ; 
      isFirst = false; 
    }
}
 
void RevMatParameter::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const 
{
  auto content = extractParameter(traln); 
  std::vector<std::string> names; 
  
  switch(content.values.size())
    {
    case 6 : 
      names = { "A<->C", "A<->G", "A<->T", "C<->G" , "C<->T", "G<->T"}; 
      break; 
    default : assert(0); 
    }

  bool isFirstG = true; 
  for(nat i = 0; i < content.values.size() ; ++i)
    {
      fileHandle << (isFirstG ? "" : "\t" ) << "r{" ;
      isFirstG = false; 
	
      bool isFirst = true; 
      for(auto &p : partitions)
	{
	  fileHandle  << (isFirst ? "": "," ) << p ; 
	  isFirst = false; 
	}	    
      fileHandle  << "}("  << names.at(i) << ")"; 
    }  
}
