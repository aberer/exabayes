#ifndef _PARTIAL_DIRICHLET_HPP 
#define _PARTIAL_DIRICHLET_HPP 

#include "proposers/AbstractProposer.hpp"
#include "model/RateHelper.hpp"


class RateDirichletProposer :  public AbstractProposer
{
public: 
  RateDirichletProposer(double minValI, double maxValI);
  
  RateDirichletProposer(const RateDirichletProposer& rhs); 

  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, log_double &hastings) ; 
  virtual AbstractProposer* clone() const ; 

private: 			// METHODS
  std::tuple<std::vector<double>, std::vector<double> > correctValues(nat whichFreq, nat numFreq, std::vector<double> newRates, std::vector<double> allNewValues); 

private: 			// ATTRIBUTES
  RateHelper _rateHelper; 
}; 

#endif
