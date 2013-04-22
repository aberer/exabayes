#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include <string>
#include "SuccessCtr.hpp" 


class AbstractProposal
{
public: 
  AbstractProposal(); 
  virtual ~AbstractProposal();

  virtual void applyToTree() = 0; 
  virtual void evaluateProposal() = 0; 
  virtual void resetTree() = 0; 
  virtual void autotune() = 0;

private: 
  string name;   
  SuccessCtr sctr; 
  // various parameters and stuff 

}; 


#endif
