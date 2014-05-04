
#include "system/BoundsChecker.hpp"
#include "TreeLengthMultiplier.hpp"
#include "math/Randomness.hpp"
#include "model/TreeAln.hpp"
#include "priors/UniformPrior.hpp"
#include "priors/AbstractPrior.hpp"

TreeLengthMultiplier::TreeLengthMultiplier( double _multiplier)
  : AbstractProposal(Category::BRANCH_LENGTHS, "TL-Mult", 1., 0.0001, 100, true)
  , multiplier(_multiplier)    
{
}


void TreeLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  storedBranches.clear(); 
  
  auto blParam = _primaryParameters[0].get(); 

  assert(_primaryParameters.size() == 1); 

  storedBranches = traln.extractBranches(blParam);

  auto  newBranches = storedBranches; 
  
  double treeScaler = rand.drawMultiplier(multiplier); 
  double initTL = 0; 
  double newTL = 0; 

  auto haveUniformPrior = dynamic_cast<UniformPrior*>(blParam->getPrior()) != nullptr;

  for(auto &b : newBranches)
    {
      auto initLength = b.getLength();
      initTL += b.getInterpretedLength(traln, blParam);

      b.setLength(  pow(initLength, treeScaler) ); 
      
      if( not BoundsChecker::checkBranch(b))
	BoundsChecker::correctBranch(b);

      double realScaling = log(b.getLength()) / log(initLength); 

      hastings *= log_double::fromAbs(realScaling); 

      double tmp = b.getInterpretedLength(traln, blParam); 
      newTL += tmp;
      
      // correct? 
      if(haveUniformPrior && blParam->getPrior()->getLogProb(ParameterContent{{tmp}}).isNegativeInfinity())
	prior.addToRatio( log_double::lowest()); 
    }

  for(auto &b : newBranches)
    traln.setBranch(b, blParam);

  if(not haveUniformPrior)
    {
      auto tmp = blParam->getPrior()->getLogProb( ParameterContent{{ newTL} }  ) / blParam->getPrior()->getLogProb( ParameterContent{{ initTL} } ); 
      prior.addToRatio( tmp ); 
    }
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

  double newParam = tuneParameter(_sctr.getBatch(), _sctr.getRatioInLastInterval(), parameter, PLL_FALSE);

  multiplier = newParam; 

  _sctr.nextBatch();
}
 
 
void TreeLengthMultiplier::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  assert(_primaryParameters.size( )== 1); 
  auto parts = _primaryParameters[0]->getPartitions(); 

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
  multiplier = cRead<decltype(multiplier)>(in); 
} 

void TreeLengthMultiplier::writeToCheckpointCore(std::ostream &out) const
{
  cWrite<decltype(multiplier)>(out, multiplier); 
} 


std::vector<nat> TreeLengthMultiplier::getInvalidatedNodes(const TreeAln& traln) const
{
  auto result = std::vector<nat>{}; 
  for(nat i = traln.getNumberOfTaxa() + 1 ; i < traln.getNumberOfNodes() + 1  ; ++i)
    result.push_back(i); 
  return result; 
} 


