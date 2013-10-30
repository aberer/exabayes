#ifndef _EXPONENTIAL_PRIOR 
#define _EXPONENTIAL_PRIOR 

#include "AbstractPrior.hpp"
#include "Density.hpp"
 
class ExponentialPrior : public AbstractPrior
{
public: 
  ExponentialPrior(double lambda) ; 
  virtual bool needsIntegration() const {return true; } 
  virtual double getLogProb(const ParameterContent& content) const  ; 
  // virtual std::vector<double> drawFromPrior(Randomness &rand)  const; 
  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const {assert(0); return ParameterContent{}; } ; 
  virtual void print(std::ostream& out ) const  ; 
  virtual double getLamda()  const  { return lambda; } 
  virtual ParameterContent getInitialValue() const; 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew )  const ; 

private: 
  double lambda; 
}; 

#endif
