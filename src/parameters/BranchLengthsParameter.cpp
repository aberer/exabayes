#include "BranchLengthsParameter.hpp"
#include "system/BoundsChecker.hpp"

void BranchLengthsParameter::applyParameter(TreeAln& traln, const ParameterContent &content) const
{
  for(auto &b : content.branchLengths)
    traln.setBranch(b, const_cast<AbstractParameter*>(dynamic_cast<const AbstractParameter* const >(this))); 
}

ParameterContent BranchLengthsParameter::extractParameter(const TreeAln &traln)  const
{
  auto result = ParameterContent{}; 
  result.branchLengths = traln.extractBranches(this); 
  return result; 
}   



void BranchLengthsParameter::verifyContent(const TreeAln &traln, const ParameterContent &content) const 
{
  for(auto &bl : content.branchLengths)
    {
      if(not BoundsChecker::checkBranch(bl)) 
	{
	  tout << "observed invalid branch " << bl << " that must not be there." << std::endl; 
	  assert(0); 
	}
    } 
} 


bool BranchLengthsParameter::priorIsFitting(const AbstractPrior &prior, const TreeAln &traln) const
{
  auto content = prior.getInitialValue();
  return content.values .size() < 2  ; 
} 

