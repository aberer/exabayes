#include "FrequencyParameter.hpp"
#include "BoundsChecker.hpp"


void FrequencyParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{  
  for(auto &m : partitions)    
    traln.setFrequencies(content.values, m); 
}


ParameterContent FrequencyParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = traln.getFrequencies(partitions[0]); 
  return result; 
}   


void FrequencyParameter::printSample(std::ostream& fileHandle, const TreeAln &traln) const 
{
  auto content =  extractParameter(traln); 
  bool isFirst = true; 
  for(auto &v : content.values)
    {
      fileHandle << (isFirst ? "" : "\t") << v ; 
      isFirst = false; 
    }
} 


void FrequencyParameter::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln)  const   
{
  auto content = extractParameter(traln); 

  std::vector<std::string> names; 
  switch(content.values.size())
    {
    case 4:       
      names = { "A" , "C", "G", "T"}; 
      break; 
    default: 
      assert(0); 
    }

  
  bool isFirstG = true; 
  for(nat i = 0; i < content.values.size() ; ++i)
    {
      fileHandle << (isFirstG ? "" : "\t" ) << "pi{" ;
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


void FrequencyParameter::verifyContent(const TreeAln &traln, const ParameterContent &content) const 
{
  auto partition = traln.getPartition(partitions[0]);
  bool ok = true; 
  ok &= BoundsChecker::checkFrequencies(content.values); 
  ok &= (content.values.size() ==  nat(partition->states)); 
  
  if(not ok)
    {
      tout << "ERROR: we had " << content << " for parameter " << this << ". Did you mis-specify a fixed prior?"  << std::endl; 
      assert(0); 
    }  
}
