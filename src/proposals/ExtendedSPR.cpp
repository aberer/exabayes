#include "ExtendedSPR.hpp"
#include "data-struct/Path.hpp"
#include "model/TreeAln.hpp"
#include "priors/AbstractPrior.hpp"
#include "TreePrinter.hpp"

// #define MOVE_INFO

// #define DEBUG_ESPR

ExtendedSPR::ExtendedSPR( double _stopProb, double _multiplier)
  : AbstractProposal(Category::TOPOLOGY,  "eSPR", 5.,  0,0, false)
  , stopProb(_stopProb) 
  , multiplier(_multiplier)    
{
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
    } while(traln.isTipNode(q) && traln.isTipNode(r) ); 
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
  boolean accepted = PLL_FALSE;   
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
		       std::make_tuple(   BranchPlain(p->number, p->back->number), 
					  BranchPlain(bla.getPrimNode(), bla.getSecNode()) ),
		       getSecondaryParameterView() ); 
  modifiedPath.clear(); 
  // END
}


/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
void ExtendedSPR::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
  auto bMode = getBranchProposalMode();
  bool multiplyBranchesUsingPosterior = bMode[0]; 
  bool outer = bMode[1]; 
  bool sequential = bMode[2]; 

  auto blParams = getBranchLengthsParameterView(); 
  drawPathForESPR( traln,rand ,stopProb); 

#ifdef MOVE_INFO
  double topoLnl = 0; 
  double oldLnl = traln.getTr()->likelihood; 
  double branchLnl = 0; 
  double impact = 0; 
#endif

  if(multiplyBranchesUsingPosterior)
    {
      auto result = move.moveBranchProposal(traln, blParams, eval, rand, outer, 0.05, sequential);
      hastings *= std::get<1>(result); 

      move.applyToTree(traln, getSecondaryParameterView() , eval, outer); 

      // BEGIN REPORT 
#ifdef MOVE_INFO
      eval.evaluate(traln, move.getEvalBranch(traln), false); 
      topoLnl = traln.getTr()->likelihood; 
#endif
      // END
      
      // save and apply the branches 
      branchesSaved = true; 
      savedBls.clear(); 
#ifdef MOVE_INFO
      impact = std::get<2>(result);
#endif
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

	      prior.addToRatio( newPr / oldPr); 
	    }

	  // set the branch
	  traln.setBranch(elem, blParams);
	}
     
#ifdef MOVE_INFO 
      // BEGIN
      for(auto n : move.getDirtyNodes(traln, outer))
	eval.markDirty(traln, n);
      eval.evaluate(traln,move.getEvalBranch(traln), false); 
      branchLnl = traln.getTrHandle().likelihood; 
      // END
#endif
    }
  else 
    {
      branchesSaved = false; 
      move.applyToTree(traln, getSecondaryParameterView() ); 
    }

#ifdef MOVE_INFO
  tout << "ESPR\t" << move.getNniDistance() 
       << "\t" << oldLnl <<  "\t" << topoLnl - oldLnl << "\t" << branchLnl - oldLnl << "\t"
       << hastings << "\t" << prior.getLnPriorRatio()<< "\t" << impact << std::endl;  
  // tout << "================================================================" << std::endl; 
#endif
}


void ExtendedSPR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion)
{  
  auto toEval = move.getEvalBranch(traln);
  auto dirtyNodes = move.getDirtyNodes(traln, false); 
  for(auto &elem : dirtyNodes)
    evaluator.markDirty( traln, elem); 


#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL " << toEval << std::endl; 
#endif
  evaluator.evaluate(traln,toEval, false); 
}


void ExtendedSPR::resetState(TreeAln &traln )
{
  auto params = getSecondaryParameterView(); 
  if(branchesSaved)
    {
      for(auto elem : savedBls)
	traln.setBranch(elem, params); 
    }

  move.revertTree(traln, getSecondaryParameterView() ); 
}


AbstractProposal* ExtendedSPR::clone() const
{
  return new ExtendedSPR( *this); 
}


std::vector<nat> ExtendedSPR::getInvalidatedNodes(const TreeAln& traln) const
{
  return move.getDirtyNodes(traln, false);
}
