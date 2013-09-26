#ifndef _EXPONENTIAL_PRIOR 
#define _EXPONENTIAL_PRIOR 

#include "AbstractPrior.hpp"
 
class ExponentialPrior : public AbstractPrior
{
public: 
  ExponentialPrior(double lambda) : lambda(lambda)
  {
  }

  virtual double getLogProb(std::vector<double> values) const 
  {
    assert(values.size() == 1); 
    double value =  values[0]; 
    double result = exponentialDensity(value, lambda); 
    result = log(result); 
    return result ; 
  }

  virtual std::vector<double> drawFromPrior(Randomness &rand)  const
  {
    double drawn = rand.drawRandExp(lambda); 
    std::vector<double> result = {drawn}; 
    return result;  
  }

  virtual void print(std::ostream& out ) const  
  {        
    out << "Exponential("  << lambda << ")" ;       
  }

  virtual double getLamda()  const  { return lambda; } 

  virtual ParameterContent getInitialValue() const
  {
    auto result = ParameterContent{}; 
    result.values.push_back( 1. / lambda); 
    return result; 
  } 


private: 
  double lambda; 
}; 

#endif
