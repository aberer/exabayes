#include "ParameterProposal.hpp"
#include "tune.h"

ParameterProposal::ParameterProposal(Category cat, std::string _name, bool modifiesBL,  
				     std::unique_ptr<AbstractProposer> _proposer, double parameter )
  : modifiesBL(modifiesBL)
  , parameter(parameter)
  , proposer(std::move(_proposer))  
{
  category = cat; 
  name = _name ; 
  relativeWeight = 0;   
} 



ParameterProposal::ParameterProposal(const ParameterProposal &rhs) 
  : AbstractProposal(rhs)  
  , modifiesBL(rhs.modifiesBL)
  , parameter(rhs.parameter)
  , proposer(std::unique_ptr<AbstractProposer>(rhs.proposer->clone()))
{
  
}


void ParameterProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  assert(primVar.size() == 1); 	// we only have one parameter to integrate over 
  // this parameter proposal works with any kind of parameters (rate
  // heterogeneity, freuqencies, revmat ... could also be extended to
  // work with AA)
  

  // extract the parameter (a handy std::vector<double> that for
  // instance contains all the frequencies)
  ParameterContent content = primVar[0]->extractParameter(traln); 
  savedContent = content; 
  
  // nasty, we have to correct for the fracchange 
  double oldFracChange = traln.getTr()->fracchange; 
  
  // we have a proposer object, that does the proposing (check out
  // ProposalFunctions.hpp) It should take care of the hastings as
  // well (to some degree )
  
  auto newValues = proposer->proposeValues(content.values, parameter, rand, hastings); 

  assert(newValues.size() == content.values.size()); 
  assert(traln.getNumBranches() == 1 ); 

  // create a copy 
  ParameterContent newContent; 
  newContent.values = newValues; 
  // use our parameter object to set the frequencies or revtmat rates
  // or what ever (for all partitions)
  primVar[0]->applyParameter(traln, newContent); 

  // now take care of the fracchange 
  double newFracChange = traln.getTr()->fracchange; 
  if(modifiesBL)
    {
      std::vector<AbstractPrior*> blPriors; 
      for(auto &v : secVar)
	blPriors.push_back(v->getPrior()); 
      prior.accountForFracChange(traln, {oldFracChange}, {newFracChange}, blPriors); 
      updateHastings(hastings, pow(newFracChange / oldFracChange, traln.getNumberOfBranches()), name); 
    }

  // a generic prior updates the prior rate 
  auto thePrior = primVar[0]->getPrior();
  prior.addToRatio(thePrior->getLogProb(newValues) - thePrior->getLogProb(savedContent.values)); 
} 


void ParameterProposal::evaluateProposal(LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior)
{
  assert(primVar.size() == 1 ); 
  evaluator->evaluatePartitions(traln, primVar[0]->getPartitions() , true); 
}
 
void ParameterProposal::resetState(TreeAln &traln, PriorBelief &prior) 
{
  assert(primVar.size() == 1 ); 
  primVar[0]->applyParameter(traln, savedContent);
  // tout << "resetting to " << savedContent << std::endl; 
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



void ParameterProposal::readFromCheckpointCore(std::istream &in)
{
  parameter = cRead<double>(in); 
} 


void ParameterProposal::writeToCheckpointCore(std::ostream &out)  const
{
  cWrite(out, parameter); 
} 
