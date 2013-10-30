#ifndef _UNIFORM_PRIOR
#define _UNIFORM_PRIOR

#include "AbstractPrior.hpp" 

class UniformPrior : public AbstractPrior
{
public: 
  UniformPrior(double minVal, double maxVal); 
  virtual double getLogProb(const ParameterContent& content)  const ; 
  // virtual std::vector<double> drawFromPrior(Randomness &rand)  const; 
  virtual bool needsIntegration() const {return true; } 
  virtual void print(std::ostream& out ) const  ; 
  virtual ParameterContent getInitialValue() const; 

  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const {assert(0); return ParameterContent{}; } ; 

  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const; 

  double getMin() const {return minVal; }
  double getMax() const {return maxVal; }

private: 
  double minVal; 
  double maxVal; 
}; 

#endif
