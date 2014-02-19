#include "BranchLengthMultiplier.hpp"
#include "TreeAln.hpp"
#include "BoundsChecker.hpp"
#include "TreeRandomizer.hpp"
#include "GibbsProposal.hpp"
#include "priors/AbstractPrior.hpp"


BranchLengthMultiplier::BranchLengthMultiplier(  double multiplier)
  : AbstractProposal(Category::BRANCH_LENGTHS, "blMult", 15., false,  0.0001,  100)
  , _multiplier(multiplier)
{
}


BranchPlain BranchLengthMultiplier::proposeBranch(const TreeAln &traln, Randomness &rand) const 
{
  if(_inSetExecution)
    return _preparedBranch;
  else  
    return determinePrimeBranch(traln,rand); 
}   


BranchPlain BranchLengthMultiplier::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
{
  return TreeRandomizer::drawBranchUniform(traln, rand); 
} 


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto b = proposeBranch(traln, rand).toBlDummy(); 

  assert(_primaryParameters.size() == 1); 
  auto param = _primaryParameters[0].get(); 

  savedBranch = traln.getBranch(b.toPlain(), param); 

  double oldZ = savedBranch.getLength();

  double
    drawnMultiplier = 0 ,
    newZ = oldZ; 

  drawnMultiplier= rand.drawMultiplier( _multiplier); 
  assert(drawnMultiplier > 0.); 
  newZ = pow( oldZ, drawnMultiplier);

  // tout << MAX_SCI_PRECISION << "proposed " << newZ <<  " from " << multiplier << " and oldBranch=" << savedBranch << std::endl; 
  
  b.setLength(newZ); 

  if(not BoundsChecker::checkBranch(b))
    BoundsChecker::correctBranch(b); 
  traln.setBranch(b, param); 

  double realMultiplier = log(b.getLength()) / log(oldZ); 
  AbstractProposal::updateHastingsLog(hastings, log(realMultiplier), _name); 

  double prNew = param->getPrior()->getLogProb( ParameterContent{{ b.getInterpretedLength(traln, param)} } ); 
  double prOld = param->getPrior()->getLogProb( ParameterContent{{ savedBranch.getInterpretedLength(traln, param) }} );
  
  prior.addToRatio(prNew - prOld );
}


void BranchLengthMultiplier::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  assert(_primaryParameters.size() == 1 ); 
  auto parts = _primaryParameters[0]->getPartitions();
  
#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL: " << savedBranch << std::endl; 
#endif
  evaluator.evaluatePartitionsWithRoot(traln,savedBranch.toPlain(), parts, false); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln) 
{
  assert(_primaryParameters.size() == 1)  ; 
  auto params = getBranchLengthsParameterView(); 
  assert(params.size() == 1); 
  auto param = params[0]; 
  traln.setBranch(savedBranch, param); 
}


void BranchLengthMultiplier::autotune() 
{
  double ratio = _sctr.getRatioInLastInterval(); 
  double newParam = tuneParameter(_sctr.getBatch(), ratio , _multiplier, false);
  _multiplier = newParam; 
  _sctr.nextBatch();
}


AbstractProposal* BranchLengthMultiplier::clone() const
{
  return new BranchLengthMultiplier(*this);
}


void BranchLengthMultiplier::readFromCheckpointCore(std::istream &in) 
{
  _multiplier = cRead<decltype(_multiplier)>(in);
} 

void BranchLengthMultiplier::writeToCheckpointCore(std::ostream &out) const
{
  cWrite<decltype(_multiplier)>(out, _multiplier); 
} 


std::pair<BranchPlain,BranchPlain> BranchLengthMultiplier::prepareForSetExecution(TreeAln &traln, Randomness &rand) 
{
  assert(_inSetExecution); 
  return std::pair<BranchPlain,BranchPlain>( 
					    determinePrimeBranch(traln,rand),
					    BranchPlain(0,0)); 
}


std::vector<nat> BranchLengthMultiplier::getInvalidatedNodes(const TreeAln &traln ) const 
{
  // auto tmp = std::vector<nat>{}; 
  // if(not traln.isTipNode(savedBranch.getPrimNode()))
  //   {
  //     auto desc1 = traln.getDescendents(savedBranch.toPlain());
  //     tmp.push_back(std::get<0>(desc1).getSecNode());
  //     tmp.push_back(std::get<1>(desc1).getSecNode());
  //   }

  // if(not traln.isTipNode(savedBranch.getSecNode()))
  //   {
  //     auto desc2 = traln.getDescendents(savedBranch.getInverted().toPlain());
  //     tmp.push_back(std::get<0>(desc2).getSecNode());
  //     tmp.push_back(std::get<1>(desc2).getSecNode());
  //   }

  // return { 
  //   savedBranch.getPrimNode() , 
  //     savedBranch.getSecNode() 
  //     } ; 
  return {}; 
}



