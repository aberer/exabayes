#include "ParameterProposal.hpp"
#include "tune.h"
#include "priors/AbstractPrior.hpp"
#include "priors/FixedPrior.hpp"
#include "BoundsChecker.hpp"

ParameterProposal::ParameterProposal(Category cat, std::string _name, bool modifiesBL,  
				     std::unique_ptr<AbstractProposer> _proposer, double parameter )
  : AbstractProposal( cat, _name)
  , modifiesBL(modifiesBL)
  , parameter(parameter)
  , proposer(std::move(_proposer))  
{
  relativeWeight = 0;   
} 



ParameterProposal::ParameterProposal(const ParameterProposal &rhs) 
  : AbstractProposal(rhs)  
  , modifiesBL(rhs.modifiesBL)
  , parameter(rhs.parameter)
  , proposer(std::unique_ptr<AbstractProposer>(rhs.proposer->clone()))
{
  
}


void ParameterProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
  auto blParams = getBranchLengthsParameterView(); 

  assert(primaryParameters.size() == 1); 	// we only have one parameter to integrate over 
  // this parameter proposal works with any kind of parameters (rate
  // heterogeneity, freuqencies, revmat ... could also be extended to
  // work with AA)
  

  // extract the parameter (a handy std::vector<double> that for
  // instance contains all the frequencies)
  auto content = primaryParameters[0]->extractParameter(traln); 
  savedContent = content; 
  
  // nasty, we have to correct for the fracchange 
  auto oldFCs = std::vector<double>{}; 
  auto newFCs = std::vector<double>{};
  for(auto &param : blParams)
    oldFCs.push_back(traln.getMeanSubstitutionRate(param->getPartitions())); 

  // we have a proposer object, that does the proposing (check out
  // ProposalFunctions.hpp) It should take care of the hastings as
  // well (to some degree)
  
  auto newValues = proposer->proposeValues(content.values, parameter, rand, hastings); 

  // tout << "proposing "  << newValues << std::endl; 
  assert(newValues.size() == content.values.size()); 

  // create a copy 
  auto newContent = ParameterContent{}; 
  newContent.values = newValues; 
  // use our parameter object to set the frequencies or revtmat rates
  // or what ever (for all partitions)
  primaryParameters[0]->applyParameter(traln, newContent); 
  
  if(modifiesBL)
    {
      for(auto &param : blParams)
	newFCs.push_back(traln.getMeanSubstitutionRate(param->getPartitions())); 

      prior.accountForFracChange(traln, oldFCs, newFCs, blParams); 

      for(nat i = 0; i < blParams.size();++i)
	{
	  auto value = traln.getNumberOfBranches() * log(newFCs.at(i) / oldFCs.at(i))   ; 
	  AbstractProposal::updateHastingsLog(hastings, value, name); 
	}
    }

  // a generic prior updates the prior rate 
  auto pr = primaryParameters[0]->getPrior();
  prior.addToRatio(pr->getLogProb(newValues) - pr->getLogProb(savedContent.values)); 
} 


void ParameterProposal::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion)
{
  auto prts = primaryParameters[0]->getPartitions(); 
#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL-CHOICE "  << branchSuggestion << std::endl; 
#endif
  evaluator.evaluatePartitionsWithRoot(traln, branchSuggestion, prts , true); 
}
 


void ParameterProposal::resetState(TreeAln &traln) 
{
  primaryParameters[0]->applyParameter(traln, savedContent);

  // for a fixed bl parameter, we have to re-scale the branch lengths after rejection again. 
  // NOTICE: this is very inefficient 
  if(modifiesBL)
    {
      for(auto &param : getBranchLengthsParameterView() )
	{
	  auto pr = dynamic_cast<FixedPrior*>(param->getPrior()); 
	  if(pr != nullptr)
	    {
	      for(auto &b : traln.extractBranches(param))
		{
		  auto content = pr->getInitialValue(); 
		  b.setConvertedInternalLength(traln,param, content.values[0]); 
		  tout << "resetting to " << b << "\t"  << b.getInterpretedLength(traln,param); 
		  if(not BoundsChecker::checkBranch(b))
		    BoundsChecker::correctBranch(b); 
		  traln.setBranch(b,param); 
		}
	    }
      	}
    }
}


void ParameterProposal::autotune()  
{
  if(not proposer->isTune() )
      return; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, not proposer->isTuneup());
  
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


std::vector<nat> ParameterProposal::getInvalidatedNodes(const TreeAln &traln) const 
{
  auto result = std::vector<nat>{}; 
  for(nat i = traln.getNumberOfTaxa() + 1 ; i < traln.getNumberOfNodes() + 1  ; ++i)
    result.push_back(i); 
  return result; 
} 
