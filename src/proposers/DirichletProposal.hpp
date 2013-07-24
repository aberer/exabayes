#ifndef DIRICHLETPROPOSAL_H
#define DIRICHLETPROPOSAL_H

#include "AbstractProposer.hpp"

class DirichletProposal : public AbstractProposer
{				
public: 
  DirichletProposal( double minVal, double maxVal); 

  DirichletProposal(const DirichletProposal& rhs) : AbstractProposer(rhs){}

  
  virtual ~DirichletProposal(){}
  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings); 
  virtual AbstractProposer* clone() const  {return new DirichletProposal(*this);  }

}; 

#endif
