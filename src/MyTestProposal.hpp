#include "AbstractProposal.hpp"




class MyTestProposal : public AbstractProposal
{
public: 
  MyTestProposal();
  virtual ~MyTestProposal();  

  virtual void applyToTree(); 
  virtual void evaluateProposal(); 
  virtual void resetTree(); 
  virtual void autotune();

private: 
  
}; 
