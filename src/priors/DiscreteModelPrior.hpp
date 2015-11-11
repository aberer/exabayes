#ifndef DISCRETE_MODEL_PRIOR
#define DISCRETE_MODEL_PRIOR

#include <unordered_map>
#include "ProtModel.hpp"
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
  virtual double getLogProb(const ParameterContent& content ) const ; 
  virtual void print(std::ostream &out) const ;

  virtual AbstractPrior* clone()const { return new  DiscreteModelPrior(*this) ; }
private: 
  std::unordered_map<ProtModel,double> _modelProbs; 
}; 

#endif
