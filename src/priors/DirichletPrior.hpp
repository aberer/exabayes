
#ifndef _DIRICHLET_PRIOR
#define _DIRICHLET_PRIOR

#include "priors/AbstractPrior.hpp"

class DirichletPrior : public AbstractPrior
{
public: 
  DirichletPrior(std::vector<double> alphas) : alphas(alphas)
  {
  }


  virtual double getLogProb(std::vector<double> values) const 
  {
    assert(values.size() == alphas.size() ); 
    double sum = 0; 
    for(auto v: values)
      sum += v; 
    assert(fabs(sum - 1.0) < 1e-6); 

    double result = densityDirichletLog(values, alphas ); 
    return result; 
  }
 
  
  virtual void print(std::ostream& out ) const 
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

  virtual std::vector<double> drawFromPrior(Randomness &rand)  const
  {
    std::vector<double> result; 
    result = rand.drawRandDirichlet(alphas); 
    return result; 
  }
  
private: 
  std::vector<double> alphas; 

} ;

#endif
