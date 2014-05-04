#include "priors/UniformPrior.hpp"


UniformPrior::UniformPrior(double minVal, double maxVal) : minVal(minVal), maxVal(maxVal)
{
}

log_double UniformPrior::getLogProb(const ParameterContent& content)  const 
{
  auto &values = content.values; 
  assert(values.size() == 1); 
  double value = values[0];

  if(minVal <= value && value <= maxVal )      
    return log_double::fromAbs(1 / (maxVal - minVal)); 
  else  
    {
      return log_double::lowest();
    }
}


// std::vector<double> UniformPrior::drawFromPrior(Randomness &rand)  const
// {
//   double val = minVal + rand.drawRandDouble01() * (maxVal - minVal); 
//   std::vector<double> result = {val}; 
//   return result; 
// }

void UniformPrior::print(std::ostream& out ) const  
{ 
  out << "Uniform("  << minVal << "," << maxVal << ")" ; 
}

ParameterContent UniformPrior::getInitialValue() const
{
  auto result = ParameterContent(); 

  assert(minVal < maxVal); 
  result.values.push_back( minVal + ( maxVal - minVal ) / 2  ) ; 
  return result; 
} 



double UniformPrior::accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double oldFc, double newFc )const 
{
  auto result = 0; 

  auto maxVal =  exp( - getMin() / newFc) , 
    minVal = exp(  - getMax() / newFc) ; 

  for(auto &b : traln.extractBranches(param))
    {
      bool less = b.getLength() < minVal; 
      bool greater = b.getLength() > maxVal; 
	      
      assert(not (less && greater)); 

      if(less || greater)
	{
	  // tout << b.getLength() << " is not okay (" << b.getInterpretedLength(traln,param) << ")" << std::endl; 
	  result += - std::numeric_limits<double>::infinity();
	}
    }
  
  return result; 
} 


double UniformPrior::getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const
{
  return 0; 
} 
