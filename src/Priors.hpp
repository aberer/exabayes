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

using namespace std; 

class AbstractPrior
{
public: 
  virtual ~AbstractPrior(){assert(0);}

  virtual vector<double> drawFromPrior(Randomness &rand)  const = 0; 
  virtual double getLogProb(vector<double> values ) const = 0; 
  virtual void print(ostream &out) const = 0;
  friend ostream& operator<<(ostream &out,  shared_ptr<AbstractPrior> rhs)
  {
    rhs->print(out); 
    return out; 
  }

}; 


 
class DirichletPrior : public AbstractPrior
{
public: 
  DirichletPrior(vector<double> alphas) : alphas(alphas)
  {
  }

  virtual double getLogProb(vector<double> values) const ; 
 
  
  virtual void print(ostream& out ) const 
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

  virtual vector<double> drawFromPrior(Randomness &rand)  const
  {
    double tmp[alphas.size()];
    for(nat i = 0; i < alphas.size(); ++i)
      tmp[i] = alphas[i]; 

    double result[alphas.size()]; 

    rand.drawRandDirichlet( result, tmp, alphas.size()); 

    vector<double> resultVect; 
    for(nat i = 0; i < alphas.size(); ++i)
      resultVect.push_back(result[i]);
    return resultVect;     
  }

private: 
  vector<double> alphas; 

} ; 


class ExponentialPrior : public AbstractPrior
{
public: 
  ExponentialPrior(double lambda) : lambda(lambda)
  {
  }

  virtual double getLogProb(vector<double> values) const 
  {
    assert(values.size() == 1); 
    double value =  values[0]; 
    double result = exponentialDensity(value, lambda); 
    result = log(result); 
    return result ; 
  }

  virtual vector<double> drawFromPrior(Randomness &rand)  const
  {
    double drawn = rand.drawRandExp(lambda); 
    vector<double> result = {drawn}; 
    return result;  
  }

  virtual void print(ostream& out ) const  
  {        
    out << "Exponential("  << lambda << ")" ;       
  }

private: 
  double lambda; 
}; 


class UniformPrior : public AbstractPrior
{
public: 
  UniformPrior(double minVal, double maxVal) : minVal(minVal), maxVal(maxVal)
  {
  }

  virtual double getLogProb(vector<double> values)  const 
  {
    assert(values.size() == 1); 
    double value = values[0];

    if(minVal < value && value < maxVal )      
      return log(1 / (maxVal - minVal)); 
    else  
      {
	double result = numeric_limits<double>::lowest(); 
	return result; 
      }
  }

  virtual vector<double> drawFromPrior(Randomness &rand)  const
  {
    double val = minVal + rand.drawRandDouble01() * (maxVal - minVal); 
    vector<double> result = {val}; 
    return result; 
  }

  virtual void print(ostream& out ) const  
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
  FixedPrior(vector<double> fixedValues) : fixedValues(fixedValues) 
  {
  } 
  
  virtual double getLogProb(vector<double> values)  const
  {    
    assert(values.size() == fixedValues.size()); 
    for(nat i = 0; i < fixedValues.size() ; ++i)
      assert(fixedValues[i] == values[i]);
    return 0; 
  }

  virtual vector<double> drawFromPrior(Randomness &rand)  const
  {
    return fixedValues; 
  }

  virtual void print(ostream &out) const 
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
  vector<double> fixedValues; 
}; 

#endif



