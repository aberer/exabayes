#ifndef _PRIORS_H
#define _PRIORS_H

#include <memory>
#include <vector>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

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

  virtual std::vector<double> drawFromPrior(Randomness &rand)  const = 0; 
  virtual double getLogProb(std::vector<double> values ) const = 0; 
  virtual void print(std::ostream &out) const = 0;
  friend std::ostream& operator<<(std::ostream &out,  AbstractPrior* rhs)
  {
    rhs->print(out); 
    return out; 
  }

}; 

 
#endif
