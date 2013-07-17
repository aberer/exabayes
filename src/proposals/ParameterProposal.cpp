#include "ParameterProposal.hpp"
#include "tune.h"

ParameterProposal::ParameterProposal(Category cat, std::string _name, bool modifiesBL,  
				     std::shared_ptr<AbstractProposer> _proposer, double parameter )
  : modifiesBL(modifiesBL)
  , parameter(parameter)
  , proposer(_proposer)  
{
  category = cat; 
  name = _name ; 
  relativeWeight = 0;   
} 


void ParameterProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  assert(primVar.size() == 1); 
  
  ParameterContent content = primVar[0]->extractParameter(traln); 
  savedContent = content; 
  
  double oldFracChange = traln.getTr()->fracchange; 
  
  auto newValues = proposer->proposeValues(content.values, parameter, rand, hastings); 
  assert(newValues.size() == content.values.size()); 
  assert(traln.getNumBranches() == 1 ); 
  
  ParameterContent newContent; 
  newContent.values = newValues; 
  primVar[0]->applyParameter(traln, newContent); 
  
  double newFracChange = traln.getTr()->fracchange; 
  if(modifiesBL)
    {
      std::vector<AbstractPrior*> blPriors; 
      for(auto &v : secVar)
	blPriors.push_back(v->getPrior()); 
      prior.accountForFracChange(traln, {oldFracChange}, {newFracChange}, blPriors); 
      updateHastings(hastings, pow(newFracChange / oldFracChange, traln.getNumberOfBranches()), name); 
    }
  
  auto thePrior = primVar[0]->getPrior();
  prior.addToRatio(thePrior->getLogProb(newValues) - thePrior->getLogProb(savedContent.values)); 
} 


void ParameterProposal::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior)
{
  assert(primVar.size() == 1 ); 
  evaluator.evaluatePartitions(traln, primVar[0]->getPartitions() ); 
}
 
void ParameterProposal::resetState(TreeAln &traln, PriorBelief &prior) 
{
  assert(primVar.size() == 1 ); 
  primVar[0]->applyParameter(traln, savedContent);
}


void ParameterProposal::autotune()  
{
  if(not proposer->isTune() )
      return; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, not proposer->isTuneup());
  
#ifdef DEBUG_PRINT_TUNE_INFO
  cout << name << ": with ratio " << sctr.getRatioInLastInterval() << ": "<< ((newParam < parameter ) ? "reducing" : "increasing") <<  "\t" << parameter << "," << newParam << endl; 
#endif
  
  parameter = newParam; 
  
  sctr.nextBatch();
}



void ParameterProposal::readFromCheckpointCore(std::ifstream &in)
{
  parameter = cRead<double>(in); 
} 


void ParameterProposal::writeToCheckpointCore(std::ofstream &out) 
{
  cWrite(out, parameter); 
} 
