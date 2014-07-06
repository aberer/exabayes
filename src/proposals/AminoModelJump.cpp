#include "AminoModelJump.hpp"
#include "BoundsChecker.hpp"
#include "Category.hpp"
#include "AbstractPrior.hpp"


AminoModelJump::AminoModelJump()
  : AbstractProposal( Category::AA_MODEL, "aaMat", 1., 0,0, true)
  , savedMod{}
  , _oldFcs{}
{
}


void AminoModelJump::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator &eval)
{
  // assert(0); 
  // must update mean substitution  rate 

  // save the old fracchanges
  auto blParams = getBranchLengthsParameterView();
  // auto oldFcs = std::vector<double>{};
  _oldFcs.clear(); 
  for(auto &blParam : blParams)
    _oldFcs.push_back(blParam->getMeanSubstitutionRate()); 

  auto primVar = getPrimaryParameterView();

  assert(primVar.size() == 1); 
  auto partitions = primVar[0]->getPartitions();

  savedMod = traln.getProteinModel(partitions[0]);

  auto newModel = savedMod; 
  while(  newModel == savedMod )
    {
      auto content = primVar[0]->getPrior()->drawFromPrior(rand ); 
      newModel = content.protModel[0];
    }

  // tout << ProtModelFun::getName(savedMod)  << " -> "<< ProtModelFun::getName(newModel) << std::endl; 

  auto newFcs = std::vector<double>{}; 
  for(auto p : partitions)
    traln.setProteinModel(p, newModel);

  // update the fracchanges
  for(auto blParam: blParams)
    blParam->updateMeanSubstRate(traln);

  for(auto blParam : blParams)
    {
      // newFcs.push_back(traln.getMeanSubstitutionRate(blParam->getPartitions()));
      newFcs.push_back(blParam->getMeanSubstitutionRate());
    }


  // prior.accountForFracChange(traln, oldFcs, newFcs, blParams); 
  auto ctr = 0u; 
  for(auto param: blParams)
    {
      prior.addToRatio( param->getPrior()->accountForMeanSubstChange(traln,param,_oldFcs.at(ctr), newFcs.at(ctr))); 
      ++ctr; 
    }
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
  auto ctr = 0; 
  for(auto &param : getBranchLengthsParameterView() )
    {
      if( not param->getPrior()->needsIntegration() )
	{
	  auto prior = param->getPrior() ;
	  for(auto &b : traln.extractBranches(param))
	    {
	      auto content = prior->getInitialValue(); 
	      b.setConvertedInternalLength(param, content.values[0]); 
	      if(not BoundsChecker::checkBranch(b))
		BoundsChecker::correctBranch(b); 
	      traln.setBranch(b,param); 
	    }
	}
      else 
	{
	  param->setMeanSubstitutionRate(_oldFcs[ctr]); 
	  ++ctr; 
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
