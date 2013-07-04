#ifndef _PRIORS_H
#define _PRIORS_H

#include <memory>
#include <vector>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "axml.h"
#include "densities.h"

#include "Randomness.hpp"

class AbstractPrior
{
public: 
  virtual ~AbstractPrior()
  {}

  virtual std::vector<double> drawFromPrior(Randomness &rand)  const = 0; 
  virtual double getLogProb(std::vector<double> values ) const = 0; 
  virtual void print(std::ostream &out) const = 0;
  friend std::ostream& operator<<(std::ostream &out,  AbstractPrior* rhs)
  {
    rhs->print(out); 
    return out; 
  }

}; 

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

private: 
  double lambda; 
}; 


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

private: 
  double minVal; 
  double maxVal; 
}; 


class FixedPrior : public AbstractPrior
{
public: 
  FixedPrior(std::vector<double> fixedValues) : fixedValues(fixedValues) 
  {
  } 
  
  virtual double getLogProb(std::vector<double> values)  const
  {    
    assert(values.size() == fixedValues.size()); 
    for(nat i = 0; i < fixedValues.size() ; ++i)
      assert(fixedValues[i] == values[i]);
    return 0; 
  }

  virtual std::vector<double> drawFromPrior(Randomness &rand)  const
  {
    return fixedValues; 
  }

  virtual void print(std::ostream &out) const 
  {
    out << "Fixed(" ;     
    bool first = true; 
    for(auto v : fixedValues)
      {
	out << (first ? "" : ",") << v ; 
	if(first) first = false; 
      }
    out << ")"; 
  }

private: 
  std::vector<double> fixedValues; 
}; 

#endif
