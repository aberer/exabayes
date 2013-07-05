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
  primVar[0]->setSavedContent(content); 
  
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
      // std::cout   << "we are modifying the bl "<< std::endl; 
      std::vector<AbstractPrior*> blPriors; 
      for(auto &v : secVar)
	blPriors.push_back(v->getPrior()); 

      // tout << std::setprecision(6) << "old frac= " << oldFracChange << "\tnew"<< newFracChange << std::endl;
      
      prior.accountForFracChange(traln, {oldFracChange}, {newFracChange}, blPriors); 
    }

  auto thePrior = primVar[0]->getPrior();
  prior.addToRatio(thePrior->getLogProb(newValues) - thePrior->getLogProb(primVar[0]->getSavedContent().values)); 
} 


void ParameterProposal::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior)
{
  assert(primVar.size() == 1 ); 
  evaluator.evaluatePartitions(traln, primVar[0]->getPartitions() ); 
}
 
void ParameterProposal::resetState(TreeAln &traln, PriorBelief &prior) 
{
  assert(primVar.size() == 1 ); 
  primVar[0]->applyParameter(traln, primVar[0]->getSavedContent());
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
