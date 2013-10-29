#ifndef _FIXE_PRIOR
#define  _FIXE_PRIOR

#include "priors/AbstractPrior.hpp"


class FixedPrior : public AbstractPrior
{
public: 
  FixedPrior(std::vector<double> fixedValues)   ; 

  virtual bool needsIntegration() const {return false; } 
  virtual double getLogProb(const ParameterContent &content )  const; 
  virtual std::vector<double> drawFromPrior(Randomness &rand )  const ; 

  virtual void print(std::ostream &out) const ; 
  virtual ParameterContent getInitialValue() const; 
 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const; 

private: 
  std::vector<double> fixedValues; 
}; 

#endif
