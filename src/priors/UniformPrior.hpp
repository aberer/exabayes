#ifndef _UNIFORM_PRIOR
#define _UNIFORM_PRIOR

#include "AbstractPrior.hpp" 

class UniformPrior : public AbstractPrior
{
public: 
  UniformPrior(double minVal, double maxVal); 
  virtual log_double getLogProb(const ParameterContent& content)  const ; 
  virtual bool needsIntegration() const {return true; } 
  virtual void print(std::ostream& out ) const  ; 
  virtual ParameterContent getInitialValue() const; 

  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const {assert(0); return ParameterContent{}; } ; 

  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const; 

  virtual AbstractPrior* clone() const { return new  UniformPrior(*this) ; }

  double getMin() const {return minVal; }
  double getMax() const {return maxVal; }

  double getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const; 

private: 
  double minVal; 
  double maxVal; 
}; 

#endif
