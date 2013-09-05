#ifndef DIRICHLETPROPOSAL_H
#define DIRICHLETPROPOSAL_H

#include "AbstractProposer.hpp"

class DirichletProposal : public AbstractProposer
{				
public: 
  /** 
      @brief constructs a dirichlet proposal
      
      @param minMaxIsRelative
      indicates whether the previous two boundary arguments are
      relative to a value (e.g., revmat) or whether they sum up to 1.
   */ 
  DirichletProposal( double minVal, double maxVal, bool minMaxIsRelative); 

  DirichletProposal(const DirichletProposal& rhs) 
    : AbstractProposer(rhs)
  {
    minMaxIsRelative = rhs.minMaxIsRelative; 
  }

  virtual ~DirichletProposal(){}
  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings); 
  virtual AbstractProposer* clone() const  {return new DirichletProposal(*this);  }

private: 
  bool minMaxIsRelative;

}; 

#endif
