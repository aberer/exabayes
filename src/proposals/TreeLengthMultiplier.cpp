#include "axml.h"


#include "BoundsChecker.hpp"
#include "TreeLengthMultiplier.hpp"
#include "Randomness.hpp"
#include "TreeAln.hpp"
#include "tune.h"
#include "priors/UniformPrior.hpp"
#include "priors/AbstractPrior.hpp"

TreeLengthMultiplier::TreeLengthMultiplier( double _multiplier)
  : AbstractProposal(Category::BRANCH_LENGTHS, "TL-Mult")
  , multiplier(_multiplier)    
{
  relativeWeight = 2 ;
  needsFullTraversal = true; 
}


void TreeLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  storedBranches.clear(); 
  
  auto blParam = primaryParameters[0].get(); 

  assert(primaryParameters.size() == 1); 

  storedBranches = traln.extractBranches(blParam);

  auto  newBranches = storedBranches; 
  
  double treeScaler = rand.drawMultiplier(multiplier); 
  double initTL = 0; 
  double newTL = 0; 
  // std::cout << "drew " <<  treeScaler << std::endl; 

  bool haveUniformPrior = dynamic_cast<UniformPrior*>(blParam->getPrior());

  for(auto &b : newBranches)
    {
      auto initLength = b.getLength();
      initTL += b.getInterpretedLength(traln, blParam);

      b.setLength(  pow(initLength, treeScaler) ); 
      
      if( not BoundsChecker::checkBranch(b))
	BoundsChecker::correctBranch(b);

      double realScaling = log(b.getLength()) / log(initLength); 
      AbstractProposal::updateHastingsLog(hastings, log(realScaling), "TL-multi"); 
      double tmp = b.getInterpretedLength(traln, blParam); 
      newTL += tmp;

      if(haveUniformPrior && blParam->getPrior()->getLogProb(ParameterContent{{tmp}})  == - std::numeric_limits<double>::infinity())
	{
	  // tout << "danger: problematic length " << tmp << std::endl; 
	  prior.addToRatio( - std::numeric_limits<double>::infinity()); 
	}
    }

  for(auto &b : newBranches)
    traln.setBranch(b, blParam);

  if(not haveUniformPrior)
    prior.addToRatio(blParam->getPrior()->getLogProb( ParameterContent{{ newTL} }  )  - blParam->getPrior()->getLogProb( ParameterContent{{ initTL} } )  ); 
}


void TreeLengthMultiplier::resetState(TreeAln &traln)  
{
  auto blParams = getPrimaryParameterView(); 
  assert(blParams.size( )== 1 ) ;
  auto blParam = blParams[0]; 
  for(auto &b : storedBranches)
    traln.setBranch(b, blParam); 
} 


void TreeLengthMultiplier::autotune()
{
  double parameter = multiplier; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, FALSE);

  multiplier = newParam; 

  sctr.nextBatch();
}
 
 
void TreeLengthMultiplier::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  assert(primaryParameters.size( )== 1); 
  auto parts = primaryParameters[0]->getPartitions(); 

#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL-CHOICE " << branchSuggestion << std::endl; 
#endif

  evaluator.evaluatePartitionsWithRoot(traln,branchSuggestion, parts, true); 
}


AbstractProposal* TreeLengthMultiplier::clone() const
{
  return new TreeLengthMultiplier(*this);
}  


void TreeLengthMultiplier::readFromCheckpointCore(std::istream &in)
{
  multiplier = cRead<double>(in); 
} 

void TreeLengthMultiplier::writeToCheckpointCore(std::ostream &out) const
{
  cWrite(out, multiplier); 
} 


std::vector<nat> TreeLengthMultiplier::getInvalidatedNodes(const TreeAln& traln) const
{
  auto result = std::vector<nat>{}; 
  for(nat i = traln.getNumberOfTaxa() + 1 ; i < traln.getNumberOfNodes() + 1  ; ++i)
    result.push_back(i); 
  return result; 
} 


