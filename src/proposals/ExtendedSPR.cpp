#include "ExtendedSPR.hpp"
#include "Path.hpp"
#include "TreeAln.hpp"
#include "priors/AbstractPrior.hpp"


#include "TreePrinter.hpp"

// #define DEBUG_ESPR

ExtendedSPR::ExtendedSPR( double _stopProb, double _multiplier)
  : AbstractProposal(Category::TOPOLOGY,  "eSPR")
  , stopProb(_stopProb) 
  , multiplier(_multiplier)    
{
  relativeWeight = 5.;
  needsFullTraversal = false; 
}

BranchPlain ExtendedSPR::determinePrimeBranch(const TreeAln &traln, Randomness &rand) const 
{
  auto  start = BranchPlain(); 
  nodeptr p,q,r; 
  do 
    {
      start = TreeRandomizer::drawBranchWithInnerNode(traln, rand); 
      
      p = start.findNodePtr(traln );
      q = p->next->back; 
      r = p->next->next->back;
    } while(traln.isTipNode(q) || traln.isTipNode(r) ); 
  return start; 
} 



// IMPORTANT TODO: is the number of branches that CAN get chosen for
// the move the same with the forward and the backward move? if not,
// this MUST be accounted for in the hastings


void ExtendedSPR::drawPathForESPR(TreeAln& traln, Randomness &rand, double stopProp )
{
  auto params =  getBranchLengthsParameterView();

  auto modifiedPath = Path{}; 
  assert(0 < stopProp && stopProp < 1.0 ); 

  assert(modifiedPath.size( ) == 0 ); 

  auto start = determinePrimeBranch(traln,rand); 
  auto p = start.findNodePtr(traln) ; 
  auto q = p->next->back; 
  auto r = p->next->next->back; 

  auto lengthR = traln.getBranch(r, params); 
  auto lengthQ = traln.getBranch(q, params); 
  traln.clipNode(q,r);
  auto newBranch = lengthR; 
  newBranch.setSecNode(q->number); 
  traln.setBranch(newBranch, params); 

  p->next->back = p->next->next->back = (nodeptr)NULL; 

  modifiedPath.append(start); 
  modifiedPath.append(BranchPlain(q->number, r->number)); 

  nodeptr currentNode = rand.drawRandDouble01() ? q : r; 
  boolean accepted = FALSE;   
  while(not accepted)
    {      
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      modifiedPath.pushToStackIfNovel(BranchPlain(currentNode->number, n->number),traln); 

      currentNode = n; 
      
      accepted = rand.drawRandDouble01() < stopProp && modifiedPath.size() > 2 ; 	
    }

  /* undo changes to the tree  */
  traln.clipNode(p->next,q); 
  newBranch = lengthQ ; 
  newBranch.setSecNode(p->number); 
  traln.setBranch(newBranch, params);

  traln.clipNode(p->next->next,r);   
  newBranch = lengthR; 
  newBranch.setSecNode(p->number); 
  traln.setBranch(newBranch, params); 
  
  /* now correct  */
  if( modifiedPath.at(2).hasNode(modifiedPath.at(1).getPrimNode()))    
    modifiedPath.at(1).setSecNode(p->number); 
  else if(modifiedPath.at(2).hasNode(modifiedPath.at(1).getSecNode()  ))    
    modifiedPath.at(1).setPrimNode( p->number); 
  else 
    assert(0); 

  /* correct the incorrectly set first branch in the path */
  modifiedPath.at(0) = modifiedPath.at(0).getThirdBranch(traln, modifiedPath.at(1)); 
  
#ifdef DEBUG_ESPR
  cout << *modifiedPath << endl; 
#endif
  
  modifiedPath.debug_assertPathExists(traln); 

  // TODO in principle, we can throw away the path of this proposal and use the move proposal  
  // move.

  // BEGIN  TODO remove modifiedPath
  auto bla = modifiedPath.at(modifiedPath.size()-1); 
  move.extractMoveInfo(traln, 
		       {  BranchPlain(p->number, p->back->number), 
			   BranchPlain(bla.getPrimNode(), bla.getSecNode()) },
		       getSecondaryParameterView() ); 
  modifiedPath.clear(); 
  // END
}


/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
void ExtendedSPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
  // auto blParam = getBranchLengthsParameterView()[0]; 
  // tout << "before move " << TreePrinter(true, true, false).printTree(traln, blParam) << std::endl;  

  // IMPORTANT 
  bool multiplyBranchesUsingPosterior = 
#ifdef PROPOSE_BRANCHES_FOR_SPR 
    true; 
#else 
  false ; 
#endif
  
  auto blParams = getBranchLengthsParameterView(); 
  // auto param = blParams[0] ; 
  // assert(blParams.size( )== 1 ); 

  drawPathForESPR( traln,rand ,stopProb); 

  // tout << "ESPR\t" << move.getNniDistance() << std::endl; 

  double hastPart = 0; 
  if(multiplyBranchesUsingPosterior)
    {  
      // NOTICE: extended! 
      
      move.proposeBranches(traln, blParams, eval, hastPart, rand, false); 
      hastings += hastPart; 
    }
  // double hastBackward =  hastPart; 

  move.applyToTree(traln, getSecondaryParameterView() ); 

  // tout << "MOVE " << move << std::endl; 
  // tout << "after move " << TreePrinter(true, true, false).printTree(traln, blParam) << std::endl;  

  // double lnlAftermove = eval.evaluate(traln,move.getEvalBranch(traln), true); 
  
  bool modifiesBl = false; 
  for(auto &v : secondaryParameters)
    modifiesBl |= v->getCategory() == Category::BRANCH_LENGTHS; 

#ifdef NO_SEC_BL_MULTI
  modifiesBl = false; 
#endif

  if(modifiesBl)
    {
      assert(0);
      // auto brPr = secondaryParameters.at(0)->getPrior();
      // move.multiplyBranches(traln, rand, hastings, prior, multiplier, {brPr} ); 
    }

  subtreeBranch = move.getSubtreeBranchAfter(traln);
  oppositeBranch = move.getOppositeBranch(traln); 
  
  // double blBefore = subtreeBranch.getInterpretedLength(traln, blParams[0]); 
  // double blAfter = 0; 

  if(multiplyBranchesUsingPosterior)
    {
      hastPart = 0; 
      auto proposedBranches = move.proposeBranches(traln, blParams, eval, hastPart, rand, true  );
      hastings += hastPart; 
      for(auto b : proposedBranches)
	{
	  // for(auto param : blParams)
	  //   {
	  //     auto bl = traln.getBranch(b.toPlain(), param); 
	  //     prior.updateBranchLengthPrior(traln, bl.getLength(), b.getLength(param), param); 
	  //   }

	  for(auto param : blParams)
	    {
	      auto dummy = traln.getBranch(b.toPlain(), param); 
	      auto prevLen = dummy.getInterpretedLength(traln, param); 
	      dummy.setLength(b.getLength(param)); 
	      auto curLen = dummy.getInterpretedLength(traln,param); 
	      
	      auto ratio = param->getPrior()->getLogProb( { curLen } ) - param->getPrior()->getLogProb( {prevLen } ); 
	      prior.addToRatio(ratio); 
	    }

	  traln.setBranch(b,blParams);
	}
    }
  // double hastForward = hastPart; 

  // tout << hastBackward << " -  " << hastForward << " => " << hastings << std::endl; 
  

  // double lnlAfterbranch = eval.evaluate(traln, move.getEvalBranch(traln), true); 
  // tout << SOME_FIXED_PRECISION; 
  // tout << "ESPR\t" << blBefore << "\t" << blAfter << "\t"  << lnlAftermove - lnlInit << "\t" << lnlAfterbranch - lnlInit << "\t" << prior.getLnPriorRatio() << "\t" << hastBackward << "\t" << hastForward << " = " << hastings <<  "\t"  << move     << std::endl; 

  // debug_checkTreeConsistency(traln); 
}


void ExtendedSPR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion)
{  
  auto toEval = move.getEvalBranch(traln);
  auto dirtyNodes = move.getDirtyNodes(); 

  // tout << "marking nodes dirty: " << dirtyNodes << std::endl; 

  for(auto &elem : dirtyNodes)
    evaluator.markDirty( traln, elem); 


#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL " << toEval << std::endl; 
#endif
  evaluator.evaluate(traln,toEval, false); 
}


void ExtendedSPR::resetState(TreeAln &traln )
{
  move.revertTree(traln, getSecondaryParameterView() ); 
  // TODO 
  auto params = getSecondaryParameterView(); 
  // auto param = params[0]; 
  // assert(params.size( )== 1 ); 
  // traln.setBranch( subtreeBranch, param);
  // traln.setBranch( oppositeBranch, param); 
}


AbstractProposal* ExtendedSPR::clone() const
{
  return new ExtendedSPR( *this); 
}


std::vector<nat> ExtendedSPR::getInvalidatedNodes(const TreeAln& traln) const
{
  return move.getDirtyNodes();
}
