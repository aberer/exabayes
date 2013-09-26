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


  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval); 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln); 
  virtual void resetState(TreeAln &traln) ; 
  virtual void autotune()  ;

  virtual AbstractProposal* clone() const {return new ParameterProposal(*this) ;   }
  virtual std::vector<nat> getInvalidatedNodes(const TreeAln &traln ) const ; 

  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  
  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0) );}

  virtual void readFromCheckpointCore(std::istream &in) ; 
  virtual void writeToCheckpointCore(std::ostream &out) const;

private: 
  bool modifiesBL; 
  double parameter; 
  std::unique_ptr<AbstractProposer> proposer;
  ParameterContent savedContent; 
}; 


#endif
