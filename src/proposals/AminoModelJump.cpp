#include "AminoModelJump.hpp"
#include "system/BoundsChecker.hpp"
#include "model/Category.hpp"
#include "priors/AbstractPrior.hpp"


AminoModelJump::AminoModelJump()
  : AbstractProposal( Category::AA_MODEL, "aaMat", 1., 0,0, true)
{
}


void AminoModelJump::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator &eval)
{
  // save the old fracchanges
  auto blParams = getBranchLengthsParameterView();
  auto oldFcs = std::vector<double>{};
  for(auto &blParam : blParams)
    oldFcs.push_back(traln.getMeanSubstitutionRate(blParam->getPartitions())); 

  auto primVar = getPrimaryParameterView();

  assert(primVar.size() == 1); 
  auto partitions = primVar[0]->getPartitions();

  savedMod = traln.getProteinModel(partitions[0]);

  auto newModel = savedMod; 
  while(  newModel == savedMod )
    {
      auto content = primVar[0]->getPrior()->drawFromPrior(rand , true ); 
      newModel = content.protModel[0];
    }

  // tout << ProtModelFun::getName(savedMod)  << " -> "<< ProtModelFun::getName(newModel) << std::endl; 

  auto newFcs = std::vector<double>{}; 
  for(auto p : partitions)
    traln.setProteinModel(p, newModel);
  for(auto blParam : blParams)
    newFcs.push_back(traln.getMeanSubstitutionRate(blParam->getPartitions()));

  prior.accountForFracChange(traln, oldFcs, newFcs, blParams); 
}


void AminoModelJump::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  evaluator.evaluate(traln,traln.getAnyBranch(), true); //TODO evaluate one partition
}

void AminoModelJump::resetState(TreeAln &traln)
{  
  auto primVar = getPrimaryParameterView();
  assert(primVar.size() == 1 ); 
  auto partitions = primVar[0]->getPartitions();
  for(auto p: partitions )
    traln.setProteinModel(p,savedMod);


  // necessary, if we have fixed branch lengths  
  for(auto &param : getBranchLengthsParameterView() )
    {
      if( not param->getPrior()->needsIntegration() )
	{
	  auto prior = param->getPrior() ;
	  for(auto &b : traln.extractBranches(param))
	    {
	      auto content = prior->getInitialValue(); 
	      b.setConvertedInternalLength(traln,param, content.values[0]); 
	      if(not BoundsChecker::checkBranch(b))
		BoundsChecker::correctBranch(b); 
	      traln.setBranch(b,param); 
	    }
	}
    }
}



void AminoModelJump::autotune() 
{
  // do nothing
}


AbstractProposal* AminoModelJump::clone() const 
{
  return new AminoModelJump(*this); 
}




std::vector<nat> AminoModelJump::getInvalidatedNodes(const TreeAln& traln) const
{
  auto result = std::vector<nat>{}; 
  for(nat i = traln.getNumberOfTaxa() + 1 ; i < traln.getNumberOfNodes() + 1  ; ++i)
    result.push_back(i); 
  return result; 
}
