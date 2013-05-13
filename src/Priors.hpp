#ifndef _PRIORS_H
#define _PRIORS_H

#include <vector>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "axml.h"
#include "densities.h"

using namespace std; 

class AbstractPrior
{
public: 
  virtual double getLogProb(vector<double> values ) const = 0; 
  virtual void print(ostream &out) const = 0;
  friend ostream& operator<<(ostream &out,  AbstractPrior *rhs)
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

  virtual double getLogProb(vector<double> values) const 
  {   
    double result = densityDirichletWrapper(values, alphas); 

    result = log(result);
    
    return result; 
  } 
  
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

private: vector<double> alphas; 


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

  virtual void print(ostream& out ) const  
  {        
    out << "Exponential("  << lambda << ")" ;       
  }

private: 
  double lambda; 
}; 





// class UniformTopology : public AbstractPrior
// {
// public: 
//   explicit UniformTopology() {}
  
//   virtual double getLogProb(vector<double>values) const
//   {
//     return 0; 
//   }

//   virtual void print(ostream &out) const
//   {
//     out << "UniformTopology()"   ; 
//   }

// } ; 


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



