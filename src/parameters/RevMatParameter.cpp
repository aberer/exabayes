#include <algorithm>
#include <functional>

#include "RevMatParameter.hpp"
#include "axml.h"

void RevMatParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  std::vector<double> tmp = content.values; 
  // tmp.insert(tmp.begin(), )
  for_each(tmp.begin(), tmp.end(), [&](double &d){return d /= *(tmp.rbegin()) ;  } );
  
  // std::cout << "setting values " ; 
  // for_each(tmp.begin(), tmp.end(), [](double d  ){std::cout << d << "," ; }) ; 
  // std::cout << std::endl; 
  
  for(auto &m : partitions)
    traln.setRevMat(tmp, m);
} 


ParameterContent RevMatParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = traln.getRevMat(partitions[0]); 
  // double sum = std::accumulate(result.values.begin(), result.values.end(), 0, std::sum¿ ); 
  double sum = 0; 
  for(auto &v : result.values )
    sum += v;  
  for_each(result.values.begin(), result.values.end(), [=] (double &d){ d/= sum; }); 
  return result; 
}   

