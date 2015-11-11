#include "GibbsBranchLength.hpp"
#include "system/BoundsChecker.hpp"
#include "GibbsProposal.hpp"
#include "priors/AbstractPrior.hpp"

GibbsBranchLength::GibbsBranchLength()
  : BranchLengthMultiplier( 0)
{
  _name = "estGibbsBL"; 
  _category = Category::BRANCH_LENGTHS; 
  _relativeWeight = 20 ; 
} 


#ifdef _EXPERIMENTAL 

#define NUM_ITER 1


#if 0 
// with two 
void GibbsBranchLength::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  // this is a bl proposal that integrates over two branch lengths. 

  // TODO i did not make sure that each pair has the same
  
  // this is still experimental 

  auto oneBranch = TreeRandomizer::drawBranchWithInnerNode(traln, rand); 

  auto desc = traln.getDescendents(oneBranch); 
  auto otherBranch = rand.drawRandDouble01() < .5 ? desc.first : desc.second;

  auto params = getBranchLengthsParameterView();
  assert(params.size() == 1 );     
  auto param = params[0]; 

  auto origBranches = std::vector<Branch>{}; 
  origBranches.push_back(traln.getBranch(oneBranch, param)); 
  origBranches.push_back(traln.getBranch(otherBranch, param)); 
  
  // update bls 
  oneBranch = origBranches[0]; 
  otherBranch = origBranches[1]; 

  // optimize both branches 
  auto oneBranchInfo = std::make_pair(0.,0.); 
  auto otherBranchInfo = std::make_pair(0.,0.); 
  for(int i = 0; i < NUM_ITER; ++i )
    {
      double f1 = 0; 
      auto optBranch = GibbsProposal::optimiseBranch(traln, oneBranch, eval, f1, oneBranchInfo.second, 30, param );
      oneBranchInfo.first = optBranch.getInterpretedLength(traln, param);
      traln.setBranch(optBranch, param); 

      optBranch = GibbsProposal::optimiseBranch(traln, otherBranch, eval, f1, otherBranchInfo.second, 30, param );
      otherBranchInfo.first = optBranch.getInterpretedLength(traln, param); 
      traln.setBranch(optBranch, param); 
    }

  auto pOneResult = GibbsProposal::propose(oneBranchInfo.first, oneBranchInfo.second, rand); 
  auto pOtherResult = GibbsProposal::propose(otherBranchInfo.first, otherBranchInfo.second, rand); 
  
  oneBranch.setConvertedInternalLength(traln, param, pOneResult[0]);
  otherBranch.setConvertedInternalLength(traln, param, pOtherResult[0]);

  if(not BoundsChecker::checkBranch(oneBranch))
    BoundsChecker::correctBranch(oneBranch, param); 
  if(not BoundsChecker::checkBranch(otherBranch))
    BoundsChecker::correctBranch(otherBranch, param); 
  
  traln.setBranch(oneBranch, param); 
  traln.setBranch(otherBranch,param); 
  
  double prevBlOne = origBranches[0].getInterpretedLength(traln,param); 
  double prevBlOther = origBranches[1].getInterpretedLength(traln,param); 

  auto hastA = logGammaDensity(prevBlOne, pOneResult[1], pOneResult[2]) - logGammaDensity(pOneResult[0], pOneResult[1], pOneResult[2]) ;  
  auto hastB = logGammaDensity(prevBlOther, pOtherResult[1], pOtherResult[2]) - logGammaDensity(pOtherResult[0], pOtherResult[1], pOtherResult[2]) ; 
  hastings += hastA ; 
  hastings += hastB; 

  prior.updateBranchLengthPrior(traln, origBranches[0].getLength(param), oneBranch.getLength(param), param);
  prior.updateBranchLengthPrior(traln, origBranches[1].getLength(param), otherBranch.getLength(param), param); 

  savedBranch = origBranches[0]; 
  extraBranch = origBranches[1]; 
}
#endif


void GibbsBranchLength::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  nat middleNode = savedBranch.getIntersectingNode(extraBranch) ;
  nat otherNode = savedBranch.getOtherNode(middleNode); 
  
  Branch b(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 

  evaluator.evaluate(traln,b, false); 
} 


void GibbsBranchLength::resetState(TreeAln &traln) 
{
  auto params = getPrimaryParameterView(); 
  assert(primaryParameters.size() == 1)  ; 
  traln.setBranch(savedBranch, params); 
  traln.setBranch(extraBranch, params);   ; 
}

#else 

void GibbsBranchLength::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto b = proposeBranch(traln, rand).toBlDummy(); 
  
  auto params = getBranchLengthsParameterView();
  assert(params.size() == 1 );     
  auto param = params[0]; 

  b = traln.getBranch(b.toPlain(),param); 
  savedBranch = b; 

  auto newBranch = GibbsProposal::drawFromEsitmatedPosterior(b, eval, traln, rand,  MAX_ITER, hastings, param); 
  traln.setBranch(newBranch, param); 

  
  double oldPr = param->getPrior()->getLogProb( ParameterContent{{ b.getInterpretedLength(traln,param) }} ); 
  double newPr = param->getPrior()->getLogProb( ParameterContent{{ newBranch.getInterpretedLength(traln,param) }} ); 

  prior.addToRatio( newPr - oldPr); 

}
#endif
