#include "StatNNI.hpp"
#include "data-struct/Path.hpp"
#include "TreeRandomizer.hpp"
#include "math/Arithmetics.hpp"
#include "AdHocIntegrator.hpp"

// #define QUICK_HACK

#if defined(QUICK_HACK ) && not defined(_EXPERIMENTAL_INTEGRATION_MODE)
#error "need integration mode "
#endif

// #define _NNI_GIBBS
// #define _EXPERIMENTAL_GIBBS
// #define _EXPERIMENTAL_GIBBS_EXTRA

#if not defined(_EXPERIMENTAL_GIBBS) &&  defined(_EXPERIMENTAL_GIBBS_EXTRA)
#error "define _EXPERIMENTAL_GIBBS as well!"
#endif


StatNNI::StatNNI( double _multiplier)
  : AbstractProposal(Category::TOPOLOGY,  "stNNI", 5., 0,0, false)
  , multiplier(_multiplier)
{
}


BranchPlain StatNNI::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const
{
  return TreeRandomizer::drawInnerBranchUniform(traln, rand); 
}


void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{    
  auto blParams = getBranchLengthsParameterView(); 
  
  auto bTmp =  determinePrimeBranch(traln, rand); 
  auto b = traln.getBranch( bTmp, blParams); 

  nodeptr p = b.findNodePtr(traln); 
  auto switchingBranch = BranchPlain( rand.drawRandDouble01() < 0.5  
				   ? p->back->next->back->number
				   : p->back->next->next->back->number, 
						      p->back->number ); 
  auto bls = std::make_tuple(b.toPlain(), switchingBranch); 
  move.extractMoveInfo(traln, bls, blParams); 

  auto priors =  std::vector<AbstractPrior* >{}; 

  auto bMode = getBranchProposalMode(); 
  auto multiplyBranches = bMode[0]; 
  auto outer = bMode[1]; 
  auto sequential = bMode[2]; 

  if(multiplyBranches )
    {
      auto result = move.moveBranchProposal(traln, blParams, eval, rand, outer, 0.05, sequential);
      hastings *= std::get<1>(result); 
      // tout << "hastings is " << result.second << std::endl;  

      move.applyToTree(traln, getSecondaryParameterView() ); 
      
      // save and apply the branches 
      branchesSaved = true; 
      savedBls.clear(); 
      for(auto elem : std::get<0>(result))
	{
	  auto curBranch = traln.getBranch(elem.toPlain(), blParams); 
	  savedBls.push_back(curBranch); 
	  
	  // inform prior 
	  for(auto param : blParams)
	    {
	      auto priorHere = param->getPrior(); 

	      auto b = elem.toBlDummy(); 
	      b.setLength(elem.getLength(param)); 
	      double lenAfterInterpret = b.getInterpretedLength(traln, param); 

	      b.setLength(curBranch.getLength(param)); 
	      double lenBeforeInterpret = b.getInterpretedLength(traln, param ); 

	      auto newPr = priorHere->getLogProb( ParameterContent{{ lenAfterInterpret}} ) ; 
	      auto oldPr = priorHere->getLogProb( ParameterContent{{  lenBeforeInterpret}}); 

	      // tout << "prRatio=" << newPr - oldPr << std::endl; 
	      // tout << "bl: " << lenBeforeInterpret << " => " << lenAfterInterpret << "\t" << newPr - oldPr << std::endl; 

	      prior.addToRatio(newPr / oldPr); 
	    }

	  // set the branch
	  traln.setBranch(elem, blParams);
	  
	}
    }
  else 
    {
      branchesSaved = false; 
      move.applyToTree(traln, getSecondaryParameterView() ); 
    }
}


void StatNNI::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion)
{
  auto evalBranch = move.getEvalBranch(traln); 

  auto dNodes = move.getDirtyNodes(traln, false); 

  // tout << SHOW(dNodes) << std::endl; 

  for (auto &elem : dNodes)
    evaluator.markDirty( traln, elem); 

#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL " << evalBranch << std::endl; 
#endif

  evaluator.evaluate(traln, evalBranch, false); 

  // tout << "END" << std::endl; 
}


void StatNNI::resetState(TreeAln &traln)  
{
  auto params = getSecondaryParameterView(); 
  if(branchesSaved)
    {
      for(auto elem : savedBls)
	traln.setBranch(elem, params); 
    }

  move.revertTree(traln, params); 
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
 
std::vector<nat> StatNNI::getInvalidatedNodes(const TreeAln& traln) const
{
  return move.getDirtyNodes(traln, false);
}
