#ifndef ALIGNMENTPROPOSAL_H
#define ALIGNMENTPROPOSAL_H

#include "AbstractProposal.hpp"
#include "proposers/AbstractProposer.hpp"
#include "LikelihoodEvaluator.hpp"

class AlignmentProposal : public AbstractProposal
{
public: 
  /////////////////
  // LIFE CYCLE  //
  /////////////////
  AlignmentProposal(Category cat, std::string name, double parameter, nat numPart, AbstractProposer* proposalPrototype, std::unique_ptr<LikelihoodEvaluator> _eval) ; 
  AlignmentProposal(const AlignmentProposal& rhs )  ; 
  AlignmentProposal& operator=(AlignmentProposal rhs); 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior); 
  virtual void autotune() ; 
  virtual AbstractProposal* clone() const ; 
  virtual void readFromCheckpointCore(std::istream &in); 
  virtual void writeToCheckpointCore(std::ostream &out) const; 

private: 
  std::vector<std::unique_ptr<AbstractProposer> > partitionProposer; 
  std::vector<double> partitionParameter; 
  std::unique_ptr<LikelihoodEvaluator> eval; 
  std::vector<ParameterContent> savedParams; 
}; 


#endif
