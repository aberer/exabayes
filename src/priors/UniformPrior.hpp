#ifndef _UNIFORM_PRIOR
#define _UNIFORM_PRIOR

#include "AbstractPrior.hpp" 

class UniformPrior : public AbstractPrior
{
public: 
  UniformPrior(double minVal, double maxVal) : minVal(minVal), maxVal(maxVal)
  {
  }

  virtual double getLogProb(std::vector<double> values)  const 
  {
    assert(values.size() == 1); 
    double value = values[0];

    if(minVal < value && value < maxVal )      
      return log(1 / (maxVal - minVal)); 
    else  
      {
	double result = std::numeric_limits<double>::lowest(); 
	return result; 
      }
  }


  virtual std::vector<double> drawFromPrior(Randomness &rand)  const
  {
    double val = minVal + rand.drawRandDouble01() * (maxVal - minVal); 
    std::vector<double> result = {val}; 
    return result; 
  }

  virtual void print(std::ostream& out ) const  
  { 
    out << "Uniform("  << minVal << "," << maxVal << ")" ; 
  }

  virtual ParameterContent getInitialValue() const
  {
    auto result = ParameterContent{}; 

    assert(minVal < maxVal); 
    result.values.push_back( minVal + ( maxVal - minVal ) / 2  ) ; 
    return result; 
  } 

private: 
  double minVal; 
  double maxVal; 
}; 

#endif
