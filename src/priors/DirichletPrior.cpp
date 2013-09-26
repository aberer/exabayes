#include "priors/DirichletPrior.hpp"



double DirichletPrior::getLogProb(std::vector<double> values) const 
{
  assert(values.size() == alphas.size() ); 
  double sum = 0; 
  for(auto v: values)
    sum += v; 
  assert(fabs(sum - 1.0) < 1e-6); 

  double result = Density::lnDirichlet(values, alphas ); 
  return result; 
}


void DirichletPrior::print(std::ostream& out ) const 
{
  out << "Dirichlet("   ; 
  bool first = true; 
  for(auto v : alphas)
    {      
      out << ( first ? "" : "," )  << v ;
      if(first)
	first = false;
    }    
  out << ")";
}



std::vector<double> DirichletPrior::drawFromPrior(Randomness &rand)  const
{
  std::vector<double> result; 
  result = rand.drawRandDirichlet(alphas); 
  return result; 
}


ParameterContent DirichletPrior::getInitialValue() const
{
  double sum = std::accumulate(alphas.begin(), alphas.end(), 0.); 
  auto result = ParameterContent{}; 
  for(auto v : alphas)
    result.values.push_back(v / sum); 
  return result; 
}
