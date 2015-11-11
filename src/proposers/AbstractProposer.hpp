#ifndef _PROPOSAL_FUNCTION_H
#define _PROPOSAL_FUNCTION_H

#include <vector>
#include <algorithm>

#include "math/Randomness.hpp"
#include "math/Density.hpp"
#include "mcmc/Chain.hpp"
#include "proposals/AbstractProposal.hpp"

//////////////
// ABSTRACT //
//////////////
class AbstractProposer
{
public:   
  AbstractProposer(bool tune, bool tuneup, double minVal, double maxVal); 
  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, log_double &hastings) = 0; 
  bool isTune() const {return _tune; } 
  bool isTuneup() const {return _tuneup; }
  void correctAbsoluteRates(std::vector<double> &values) const ; 
  virtual AbstractProposer* clone() const = 0; 

protected: 
  bool _tune; 
  bool _tuneup; 
  double _minVal; 
  double _maxVal; 
}; 


#endif


