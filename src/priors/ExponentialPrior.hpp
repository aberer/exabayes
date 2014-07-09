#ifndef _EXPONENTIAL_PRIOR 
#define _EXPONENTIAL_PRIOR 

#include "AbstractPrior.hpp"
#include "math/Density.hpp"
 
class ExponentialPrior : public AbstractPrior
{
public: 
  ExponentialPrior(double lambda) ; 
  virtual bool needsIntegration() const {return true; } 
  virtual log_double getLogProb(const ParameterContent& content) const  ; 
  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const {assert(0); return ParameterContent{}; } ; 
  virtual void print(std::ostream& out ) const  ; 
  virtual double getLamda()  const  { return lambda; } 
  virtual ParameterContent getInitialValue() const; 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew )  const ; 

  virtual double getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const; 

  virtual AbstractPrior* clone() const { return new ExponentialPrior(*this) ; }

private: 
  double lambda; 
}; 

#endif
