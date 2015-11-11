#ifndef _PRIORS_H
#define _PRIORS_H

#include <vector>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>

#include "math/Density.hpp"

#include "math/Randomness.hpp"
#include "model/TreeAln.hpp"

#include "parameters/ParameterContent.hpp"


class AbstractPrior
{
public: 
  AbstractPrior()
    : _keepInitData{false}
  { }
  
  virtual ~AbstractPrior() {}

  /** 
      @brief obtains a pre-defined initial value, depending on the
      prior.
  */
  virtual ParameterContent getInitialValue() const = 0; 
  virtual double accountForMeanSubstChange( TreeAln &traln, const AbstractParameter* param , double myOld, double myNew ) const = 0; 

  virtual ParameterContent drawFromPrior(Randomness &rand, bool uniform)  const = 0 ; 
  virtual log_double getLogProb( const ParameterContent &content ) const = 0; 
  virtual void print(std::ostream &out) const = 0;

  virtual bool needsIntegration() const = 0; 

  virtual AbstractPrior* clone() const = 0;

  virtual double getFirstDerivative(const TreeAln &traln, const AbstractParameter& param) const = 0; 

  friend std::ostream& operator<<(std::ostream &out,  AbstractPrior* rhs)
  {
    rhs->print(out); 
    return out; 
  }
  
  bool isKeepInitData() const {return _keepInitData; }
  void setKeepInitData()  { _keepInitData = true; }

protected: 
  bool _keepInitData;

}; 

 
#endif
