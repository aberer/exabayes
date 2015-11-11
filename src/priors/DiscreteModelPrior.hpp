#ifndef DISCRETE_MODEL_PRIOR
#define DISCRETE_MODEL_PRIOR

#include <unordered_map>
#include "model/ProtModel.hpp"
#include "priors/AbstractPrior.hpp"

class DiscreteModelPrior : public AbstractPrior
{
public: 
  DiscreteModelPrior(std::unordered_map<ProtModel,double> model); 
  
  // if we have only one model this is basically a fixed prior 
  virtual bool needsIntegration() const {assert(_modelProbs.size() > 0 ); return _modelProbs.size() > 1 ; } 

  virtual ParameterContent getInitialValue() const ; 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param , double myOld, double myNew ) const ; 
  ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const ; 
  virtual log_double getLogProb(const ParameterContent& content ) const ; 
  virtual void print(std::ostream &out) const ;

  double getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const {assert(0); return 0; } // doesnt have that 

  virtual AbstractPrior* clone()const { return new  DiscreteModelPrior(*this) ; }
private: 
  std::unordered_map<ProtModel,double> _modelProbs; 
}; 

#endif
