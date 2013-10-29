#include "priors/ExponentialPrior.hpp"

ExponentialPrior::ExponentialPrior(double lambda) : lambda(lambda)
{
}

double ExponentialPrior::getLogProb(const ParameterContent& content) const 
{
  auto& values = content.values; 
  assert(values.size() == 1); 
  return Density::lnExponential(values[0], lambda); 
}

std::vector<double> ExponentialPrior::drawFromPrior(Randomness &rand)  const
{
  double drawn = rand.drawRandExp(lambda); 
  std::vector<double> result = {drawn}; 
  return result;  
}

void ExponentialPrior::print(std::ostream& out ) const  
{        
  out << "Exponential("  << lambda << ")" ;       
}



double ExponentialPrior::accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const 
{
  double blInfluence = 0; 
  for(auto &b : traln.extractBranches(param))
    blInfluence += log(b.getLength());

  double result = (myNew - myOld) * lambda * blInfluence; 

  return result; 
}



ParameterContent ExponentialPrior::getInitialValue() const
{
  auto result = ParameterContent{}; 
  result.values.push_back( 1. / lambda); 
  return result; 
} 
