#include "AminoModelJump.hpp"
#include "Category.hpp"
#include "priors/AbstractPrior.hpp"


AminoModelJump::AminoModelJump()
  : AbstractProposal( Category::AA_MODEL, "aaMat")
{
  relativeWeight = 1. ;
}


void AminoModelJump::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator &eval)
{
  // save the old fracchanges
  auto blParams = getBranchLengthsParameterView();
  auto oldFcs = std::vector<double>{};
  for(auto &blParam : blParams)
    oldFcs.push_back(traln.getMeanSubstitutionRate(blParam->getPartitions())); 

  auto primVar = getPrimaryParameterView();

  assert(primVar.size() == 1); 
  auto partitions = primVar[0]->getPartitions();

  savedMod = traln.getModelAssignment(partitions[0]);

  auto newModel = savedMod; 
  while(  newModel == savedMod )
    {
      auto content = primVar[0]->getPrior()->drawFromPrior(rand , true ); 
      newModel = content.protModel[0];
    }

  // tout << "new model is "<< ProtModelFun::getName(newModel)<< " old was "<< ProtModelFun::getName(savedMod) << std::endl; 
  // auto freqs =  traln.getFrequencies(partitions[0]);
  // tout << "frequencies are " << freqs << std::endl;  

  auto newFcs = std::vector<double>{}; 
  for(auto p : partitions)
    traln.setModelAssignment(p, newModel);
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
    traln.setModelAssignment(p,savedMod);
}



void AminoModelJump::autotune() 
{
  // do nothing
}


AbstractProposal* AminoModelJump::clone() const 
{
  return new AminoModelJump(*this); 
}

