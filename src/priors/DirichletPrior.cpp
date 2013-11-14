#include "priors/DirichletPrior.hpp"
#include <numeric>


double DirichletPrior::getLogProb(const ParameterContent& content) const 
{
  auto &values = content.values; 
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


ParameterContent DirichletPrior::getInitialValue() const
{
  double sum = std::accumulate(alphas.begin(), alphas.end(), 0.); 
  auto result = ParameterContent{}; 
  for(auto v : alphas)
    result.values.push_back(v / sum); 
  return result; 
}
