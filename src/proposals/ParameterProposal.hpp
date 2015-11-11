#ifndef _PARAMETER_PROPOSALS 
#define _PARAMETER_PROPOSALS 

#include <memory> 

#include "TreeAln.hpp"
#include "ProposalFunctions.hpp"
#include "AbstractProposal.hpp"
#include "ParameterProposal.hpp"


class ParameterProposal : public AbstractProposal
{
public: 
  ParameterProposal(Category cat, std::string _name, bool modifiesBL, std::shared_ptr<AbstractProposer> _proposer, double parameter ); 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand); 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior); 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune()  ;

  virtual AbstractProposal* clone() const {return new ParameterProposal(*this) ;   }

private: 
  bool modifiesBL; 
  double parameter; 
  std::shared_ptr<AbstractProposer> proposer; 
  ParameterContent savedContent; 
}; 


#endif
