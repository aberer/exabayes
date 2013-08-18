#ifndef _BRANCH_INTEGRATOR
#define _BRANCH_INTEGRATOR

#include "BranchLengthMultiplier.hpp"

class BranchIntegrator  : public BranchLengthMultiplier
{
public: 
  BranchIntegrator(double _mult )
    : BranchLengthMultiplier(_mult)
  {
    name =  "blInt"; 
    relativeWeight = 20; 
    this->category = Category::BRANCH_LENGTHS;     
  }

  virtual Branch proposeBranch(const TreeAln &traln, Randomness &rand) const
  {
    
    return toPropose; 
  } 
  
  void setToPropose(Branch b) { toPropose = b; }

  virtual AbstractProposal* clone() const {return new BranchIntegrator(*this); }   

  virtual void autotune(){}

private:
  Branch toPropose; 
  

}; 


#endif
