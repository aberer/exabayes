#include <string>


class AbstractProposal
{
public: 
  virtual ~AbstractProposal();

  virtual void applyToTree(); 
  virtual void evaluateProposal(); 
  virtual void resetTree(); 
  virtual void autotune();

  
private: 
  string name; 
  // various parameters and stuff 

}; 
