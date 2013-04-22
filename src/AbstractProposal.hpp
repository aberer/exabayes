#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include <string>
#include "SuccessCtr.hpp" 


class AbstractProposal
{
public: 
  AbstractProposal(double weight); 
  virtual ~AbstractProposal();

  virtual void applyToTree(); 
  virtual void evaluateProposal(); 
  virtual void resetTree(); 
  virtual void autotune();

private: 
  string name; 
  
  // various parameters and stuff 

}; 


#endif
