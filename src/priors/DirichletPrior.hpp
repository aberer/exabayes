
#ifndef _DIRICHLET_PRIOR
#define _DIRICHLET_PRIOR

#include "priors/AbstractPrior.hpp"

class DirichletPrior : public AbstractPrior
{
public: 
  DirichletPrior(std::vector<double> alphas) : alphas(alphas)
  {
  }

  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const {assert(0); return ParameterContent{}; } ; 

  virtual log_double getLogProb( const ParameterContent& content) const ; 
  virtual void print(std::ostream& out ) const ; 

  virtual ParameterContent getInitialValue() const; 

  virtual bool needsIntegration() const {return true; } 

  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const {return 0 ; } 

  virtual AbstractPrior* clone() const { return new  DirichletPrior(*this) ; }

  double getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const {assert(0); return 0; } // doesnt have that
  
private: 
  std::vector<double> alphas; 

} ;

#endif
