#ifndef _PARAMETER_PROPOSALS 
#define _PARAMETER_PROPOSALS 

#include <memory> 

#include "TreeAln.hpp"
#include "AbstractProposal.hpp"
#include "proposers/AbstractProposer.hpp"
#include "ParameterProposal.hpp"



class ParameterProposal : public AbstractProposal
{
public: 
  ParameterProposal(Category cat, std::string _name, bool modifiesBL, std::unique_ptr<AbstractProposer> _proposer, double parameter ); 
  ParameterProposal(const ParameterProposal &prop); 


  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand); 
  virtual void evaluateProposal(LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior); 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune()  ;

  virtual AbstractProposal* clone() const {return new ParameterProposal(*this) ;   }

  virtual void readFromCheckpointCore(std::istream &in) ; //  { in >> parameter; readDelimiter(in);    } 
  virtual void writeToCheckpointCore(std::ostream &out) const; //   {out << parameter << DELIM; } 

private: 
  bool modifiesBL; 
  double parameter; 
  std::unique_ptr<AbstractProposer> proposer;
  ParameterContent savedContent; 
}; 


#endif
