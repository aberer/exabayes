#ifndef MULTIPLIERPROPOSAL_H
#define MULTIPLIERPROPOSAL_H

#include "AbstractProposer.hpp"

////////////////
// MULTIPLIER //
////////////////
class MultiplierProposal : public AbstractProposer
{
public: 
  MultiplierProposal(double minVal, double maxVal); 
  
  MultiplierProposal(const MultiplierProposal& rhs): AbstractProposer(rhs) {}

  virtual ~MultiplierProposal() {}

  virtual AbstractProposer* clone() const  {return new MultiplierProposal(*this);  }


  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings); 
}; 

#endif
