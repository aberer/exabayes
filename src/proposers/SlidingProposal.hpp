#ifndef SLIDINGPROPOSAL
#define SLIDINGPROPOSAL

#include "AbstractProposer.hpp"


class SlidingProposal : public AbstractProposer
{
public: 
  SlidingProposal(double minVal, double maxVal); 
  virtual ~SlidingProposal(){}

  SlidingProposal(const SlidingProposal &rhs) : AbstractProposer(rhs){}

  double proposeOneValue(double oldVal, double parameter, Randomness &rand, double &hastings); 

  std::vector<double> proposeRelativeMany(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings); 


  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings); 

  virtual AbstractProposer* clone() const  {return new SlidingProposal(*this);  }

}; 

#endif
