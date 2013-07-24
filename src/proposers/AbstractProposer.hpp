#ifndef _PROPOSAL_FUNCTION_H
#define _PROPOSAL_FUNCTION_H

#include <vector>
#include <algorithm>

#include "Randomness.hpp"
#include "densities.h"
#include "Chain.hpp"
#include "AbstractProposal.hpp"

//////////////
// ABSTRACT //
//////////////
class AbstractProposer
{
public:   
  virtual ~AbstractProposer(){}

  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings) = 0; 

  bool isTune() const {return tune; } 
  bool isTuneup() const {return tuneup; }

  virtual AbstractProposer* clone() const = 0; 
  // AbstractProposer* clone() const {return new AbstractProposer(*this) ; } 

  // TODO would be cool 
  // void setTunedParameter(double param ) { tunedParameter = param; }

protected: 
  bool tune; 
  bool tuneup; 
  double minVal; 
  double maxVal; 
}; 


#endif



// NOT USE YET 


// disabling those for now. Just more stuff to maintain. 
// we can easily reactivate them 
#if 0 
class ExponentialProposal : public AbstractProposer
{
public: 
  ExponentialProposal()
  {
    tune = false; 
    tuneup = false; 
  }

  virtual std::vector<double> proposeValues(std::vector<double> oldValue, double parameter, Randomness &rand, double &hastings)
  {
    // TODO @kassian: how to modify the hastings? 
    // assert(0); 
    return std::vector<double>();

  }

};


class BiunifProposal : public AbstractProposer
{  
  // TODO incorporate the parametr 
public: 
  BiunifProposal()
  {
    tune = true; 
    tuneup = true; 
  }

  virtual std::vector<double> proposeValues(std::vector<double> oldValue, double parameter, Randomness &rand, double &hastings)
  {
    // TODO @ kassian: how to modify the hastings?  
    assert(0); 
    return std::vector<double>();
  }

} ;  

#endif
