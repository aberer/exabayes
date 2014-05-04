#include "priors/ExponentialPrior.hpp"
#include "parameters/AbstractParameter.hpp"
#include "parameters/BranchLengthsParameter.hpp"

ExponentialPrior::ExponentialPrior(double lambda) : lambda(lambda)
{
}

log_double ExponentialPrior::getLogProb(const ParameterContent& content) const 
{
  auto& values = content.values; 
  assert(values.size() == 1); 
  return Density::lnExponential(values[0], lambda); 
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


double ExponentialPrior::getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const
{
  // this assert is clearly ridiculous; however it is not worth making
  // the situation even more complicated...
  assert( dynamic_cast<BranchLengthsParameter*>(const_cast<AbstractParameter*>(&param)) != nullptr ); 
  auto fracchange = traln.getMeanSubstitutionRate(param.getPartitions()); 

  // return 1; 
  return fracchange * lambda; 
}
