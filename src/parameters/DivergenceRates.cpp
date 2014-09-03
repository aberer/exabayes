#include "DivergenceRates.hpp"

#include "Category.hpp"


DivergenceRates::DivergenceRates(nat id, nat idOfMyKind, std::vector<nat> partitions, nat numberOfTaxa)
  : AbstractParameter(Category::DIVERGENCE_RATES, id, idOfMyKind, partitions, 0)
  , _rateAssignments(numberOfTaxa-1, 0 )
  , _rates(1,1.)
{
}



void DivergenceRates::applyParameter(TreeAln& traln,  const ParameterContent &content) const 
{
}


ParameterContent DivergenceRates::extractParameter(const TreeAln &traln)  const  
{
  assert(0); 
  return ParameterContent{{}}; 
}

   
void DivergenceRates::printSample(std::ostream& fileHandle, const TreeAln &traln ) const 
{
}
 
void DivergenceRates::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  
{
}
 
void DivergenceRates::verifyContent(const TreeAln &traln, const ParameterContent &content) const 
{
  
}


log_double DivergenceRates::getPriorValue(const TreeAln& traln) const
{
  auto result = log_double::fromAbs(1); 
  
  // that only partly makes sense ... 
  // result *= _prior->getLogProb( _rates);

  return result; 
} 



void DivergenceRates::setPrior(const std::unique_ptr<AbstractPrior> &prior) 
{
  // meh meh meh 

  AbstractParameter::setPrior(prior); 
  _rates.resize(prior->getInitialValue().values.size());
}
