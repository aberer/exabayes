#include <algorithm>
#include <functional>

#include "comm/ParallelSetup.hpp"

#include "RateHelper.hpp"
#include "DnaAlphabet.hpp"
#include "AminoAcidAlphabet.hpp"

#include "BoundsChecker.hpp"
#include "GlobalVariables.hpp"

#include "RevMatParameter.hpp"
#include "axml.h"

void RevMatParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  auto tmp = content.values; 
  RateHelper::convertRelativeToLast(tmp); 
  
  for(auto &m : _partitions)
    traln.setRevMat(tmp, m);

  // auto newStuff = traln.getRevMat(_partitions[0]);
  // tout << "applying new values " ;  
  // std::copy(newStuff.begin(), newStuff.end(), std::ostream_iterator<double>(tout, ",")); 
  // tout << std::endl;   
  // tout << "fracchange = " << traln.getTr()->fracchange << std::endl; 
} 


ParameterContent RevMatParameter::extractParameter(const TreeAln &traln )  const
{
  ParameterContent result; 
  result.values = traln.getRevMat(_partitions.at(0)); 
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
  
  auto size = content.values.size(); 
  switch(size)
    {
    case 6 : 
      names = DnaAlphabet().getCombinations();
      break; 
    case 190: 			// really? 
      names = AminoAcidAlphabet().getCombinations(); 
      break; 
    default : 
      {
	std::cerr << "error in revMatParameter: encountered case >" <<  size << std::endl; 
	assert(0); 
      }
    }

  bool isFirstG = true; 
  for(nat i = 0; i < content.values.size() ; ++i)
    {
      fileHandle << (isFirstG ? "" : "\t" ) << "r{" ;
      isFirstG = false; 
	
      bool isFirst = true; 
      for(auto &p : _partitions)
	{
	  fileHandle  << (isFirst ? "": "," ) << p ; 
	  isFirst = false; 
	}	    
      fileHandle  << "}("  << names.at(i) << ")"; 
    }  
}


void RevMatParameter::verifyContent(const TreeAln&traln,  const ParameterContent &content) const 
{
  auto& partition = traln.getPartition(_partitions[0]); 
  auto num = RateHelper::numStateToNumInTriangleMatrix(partition.states);

  bool ok = true; 

  ok &= content.values.size( )== num ; 
  
  auto newValues = content.values; 
  RateHelper::convertRelativeToLast(newValues); 

  // for(auto &v : newValues)
  //   v /= *(newValues.rbegin()); 
  
  if(ok)
    ok &= BoundsChecker::checkRevmat(newValues); 

  if(not ok )
    {
      tout << "Wrong content " << content << " for parameter "
      << this << ". Did you mis-specify a fixed prior or are your input values to extreme?" << std::endl; 
      assert(0); 
    }
} 


void RevMatParameter::checkSanityPartitionsAndPrior(const TreeAln &traln) const 
{
  checkSanityPartitionsAndPrior_FreqRevMat(traln);
  nat numStates = traln.getPartition(_partitions[0]).states; 
  
  auto initVal = _prior->getInitialValue(); 
  
  if( initVal.values.size() != RateHelper::numStateToNumInTriangleMatrix(numStates) && 
      initVal.protModel.size() == 0 )
    {
      tout << "Error while processing parsed priors: you specified prior " << _prior.get() << " for parameter "; 
      printShort(tout) << " that is not applicable." << std::endl; 
      ParallelSetup::genericExit(-1); 
    }
}
