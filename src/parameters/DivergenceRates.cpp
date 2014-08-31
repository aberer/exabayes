#include "DivergenceRates.hpp"

#include "Category.hpp"


DivergenceRates::DivergenceRates(nat id, nat idOfMyKind, std::vector<nat> partitions, nat numberOfTaxa, nat numRateCats)
  : AbstractParameter(Category::DIVERGENCE_TIMES, id, idOfMyKind, partitions, 0)
  , _rateAssignments(numberOfTaxa-1, 0 )
  , _rates(numRateCats,1.)
{
}



void DivergenceRates::applyParameter(TreeAln& traln,  const ParameterContent &content) const 
{
}

 
ParameterContent DivergenceRates::extractParameter(const TreeAln &traln)  const  
{
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
 
