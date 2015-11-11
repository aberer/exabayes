#ifndef _FIXE_PRIOR
#define  _FIXE_PRIOR

#include "priors/AbstractPrior.hpp"


class FixedPrior : public AbstractPrior
{
public: 
  FixedPrior(std::vector<double> fixedValues)   ; 

  virtual bool needsIntegration() const {return false; } 
  virtual log_double getLogProb(const ParameterContent &content )  const; 

  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const {assert(0); return ParameterContent{}; } ; 

  virtual void print(std::ostream &out) const ; 
  virtual ParameterContent getInitialValue() const; 

  virtual AbstractPrior* clone() const { return new  FixedPrior(*this) ; }
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const; 

  double getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const {assert(0); return 0; }

private: 
  std::vector<double> fixedValues; 
}; 

#endif
