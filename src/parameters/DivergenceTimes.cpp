#include "DivergenceTimes.hpp"
#include "Category.hpp"

DivergenceTimes::DivergenceTimes(nat id, nat idOfMyKind, std::vector<nat> partitions, NodeAge age)
  : AbstractParameter(Category::DIVERGENCE_TIMES, id, idOfMyKind, partitions, 0 )
  , _nodeAge{age}
{
}


void DivergenceTimes::applyParameter(TreeAln& traln,  const ParameterContent &content)
{
}

 
ParameterContent DivergenceTimes::extractParameter(const TreeAln &traln)  const  
{
  assert(0); 
  return ParameterContent{{}} ; 
}


void DivergenceTimes::printSample(std::ostream& fileHandle, const TreeAln &traln ) const 
{
}

 
void DivergenceTimes::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  
{
}

 
void DivergenceTimes::verifyContent(const TreeAln &traln, const ParameterContent &content) const
{
} 





log_double DivergenceTimes::getPriorValue(const TreeAln& traln) const
{
  // assert(0); return log_double::fromAbs(1); 
  // TODO should extract all branches and evaluate the prior ... not doing this here. .. 
  return log_double::fromAbs(1.);
} 
