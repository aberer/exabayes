#include <algorithm>
#include <functional>


#include "BoundsChecker.hpp"
#include "GlobalVariables.hpp"

#include "RevMatParameter.hpp"
#include "axml.h"

void RevMatParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  auto tmp = content.values; 

  for_each(tmp.begin(), tmp.end(), [&](double &d){return d /= *(tmp.rbegin()) ;  } );
  
  for(auto &m : partitions)
    traln.setRevMat(tmp, m);

  auto newStuff = traln.getRevMat(partitions[0]);
  
  // tout << "applying new values " ;  
  // std::copy(newStuff.begin(), newStuff.end(), std::ostream_iterator<double>(tout, ",")); 
  // tout << std::endl;   
  // tout << "fracchange = " << traln.getTr()->fracchange << std::endl; 
} 


ParameterContent RevMatParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = traln.getRevMat(partitions.at(0)); 
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





void RevMatParameter::verifyContent(const TreeAln&traln, const ParameterContent &content) const 
{
  auto partition = traln.getPartition(partitions[0]); 
  auto num = numStateToNumInTriangleMatrix(partition->states);
  
  // auto sum = std::accumulate(content.values.begin(), content.values.end(), 0.); 

  bool ok = true; 
  // ok &= fabs(sum - 1.0 )  < 1e-2 ; 

  ok &= content.values.size( )== num ; 

  auto newValues = content.values; 
  for(auto &v : newValues)
    v /= *(newValues.rbegin()); 
  
  if(ok)
    ok &= BoundsChecker::checkRevmat(newValues); 

  if(not ok )
    {
      tout << "Wrong content " << content << " for parameter "
      << this << ". Did you mis-specify a fixed prior or are your input values to extreme?" << std::endl; 
      assert(0); 
    }
} 






