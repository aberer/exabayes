
#ifndef _DIRICHLET_PRIOR
#define _DIRICHLET_PRIOR

#include "priors/AbstractPrior.hpp"

class DirichletPrior : public AbstractPrior
{
public: 
  DirichletPrior(std::vector<double> alphas) : alphas(alphas)
  {
  }

  virtual double getLogProb(std::vector<double> values) const ; 
  virtual void print(std::ostream& out ) const ; 
  virtual std::vector<double> drawFromPrior(Randomness &rand)  const; 
  virtual ParameterContent getInitialValue() const; 

  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param, double myOld, double myNew ) const {return 0 ; } 
  
private: 
  std::vector<double> alphas; 

} ;

#endif
