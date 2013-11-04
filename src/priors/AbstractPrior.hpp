#ifndef _PRIORS_H
#define _PRIORS_H

#include <vector>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>

#include "axml.h"
#include "Density.hpp"

#include "Randomness.hpp"

#include "parameters/ParameterContent.hpp"

class AbstractPrior
{
public: 
  virtual ~AbstractPrior()
  {}

  /** 
      @brief obtains a pre-defined initial value, depending on the
      prior.
  */
  virtual ParameterContent getInitialValue() const = 0; 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param , double myOld, double myNew ) const = 0; 

  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const = 0 ; 
  virtual double getLogProb( const ParameterContent &content ) const = 0; 
  virtual void print(std::ostream &out) const = 0;

  virtual bool needsIntegration() const = 0; 

  virtual AbstractPrior* clone() const = 0;

  friend std::ostream& operator<<(std::ostream &out,  AbstractPrior* rhs)
  {
    rhs->print(out); 
    return out; 
  }

}; 

 
#endif
