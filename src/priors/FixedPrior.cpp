#include "priors/FixedPrior.hpp"
#include "parameters/AbstractParameter.hpp"
#include "Category.hpp"
#include "BoundsChecker.hpp"

FixedPrior::FixedPrior(std::vector<double> _fixedValues)  : fixedValues(_fixedValues) 
{
}


double FixedPrior::getLogProb(const ParameterContent& content)  const
{    
  auto &values = content.values; 

  // branch lengths case! 
  if(fixedValues.size( )== 1 )
    {
      if(fabs (values[0] - fixedValues[0] ) > 1e-6)
	{
	  tout << "Fatal: for fixed prior, value should be "   << fixedValues[0] << " but actually is " << values[0] << std::endl; 
	  assert(0); 
	}
    }
  else 
    {
      for(nat i = 0; i < fixedValues.size() ; ++i)
	assert(fixedValues[i] == values[i]);
    }

  return 0; 
}


ParameterContent FixedPrior::getInitialValue() const
{
  auto result = ParameterContent{}; 
  result.values = fixedValues; 
  return result; 
} 


std::vector<double> FixedPrior::drawFromPrior(Randomness &rand)  const 
{
  return fixedValues; 
}


void FixedPrior::print(std::ostream &out) const 
{
  out << "Fixed(" ;     
  bool first = true; 
  for(auto v : fixedValues)
    {
      out << (first ? "" : ",") << v ; 
      if(first) first = false; 
    }
  out << ")"; 
}


double FixedPrior::accountForMeanSubstChange(TreeAln  &traln, const AbstractParameter* param, double myOld, double myNew ) const 
{
  if(param->getCategory() == Category::BRANCH_LENGTHS)
    {
      // first check, if re-scaling is possible 
      auto bls = traln.extractBranches(param); 

      for(auto &b : bls)
	b.setConvertedInternalLength(traln, param, fixedValues[0]);
      
      if( std::any_of(bls.begin(), bls.end(), 
		      [](BranchLength &b){ return not BoundsChecker::checkBranch(b) ;  }) ) 
	return - std::numeric_limits<double>::infinity(); 
      
      for(auto &b : bls)
	traln.setBranch(b,param); 

      return 0; 
    }
  else 
    {
      assert(0); 
      return 0; 
    }
}
