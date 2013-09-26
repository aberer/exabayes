#ifndef _FIXE_PRIOR
#define  _FIXE_PRIOR

#include "AbstractPrior.hpp"


class FixedPrior : public AbstractPrior
{
public: 
  FixedPrior(std::vector<double> fixedValues)   ; 

  virtual double getLogProb(std::vector<double> values)  const; 
  virtual std::vector<double> drawFromPrior(Randomness &rand)  const ; 

  virtual void print(std::ostream &out) const ; 
  virtual ParameterContent getInitialValue() const; 
 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const; 

private: 
  std::vector<double> fixedValues; 
}; 

#endif
