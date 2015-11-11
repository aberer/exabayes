#include "ExtendedTBR.hpp"
#include "model/Category.hpp"
#include "TreeRandomizer.hpp"


// TODO the disorient is still  very inefficient 

ExtendedTBR::ExtendedTBR( double extensionProb, double multiplier)
  : AbstractProposal(Category::TOPOLOGY , "eTBR", 5., 0,0, false )
  , _extensionProbability(extensionProb)
  , _multiplier(multiplier)
{
}


void ExtendedTBR::buildPath(Path &path, BranchPlain bisectedBranch, TreeAln &traln, Randomness &rand )
{
  auto params = getBranchLengthsParameterView();

  double stopProb = _extensionProbability; 

  nodeptr p= bisectedBranch.findNodePtr(traln ); 
  path.clear();

  nodeptr pn = p->next->back,
    pnn = p->next->next->back; 
  
  
  auto pnBranch = traln.getBranch(pn, params),
    pnnBranch = traln.getBranch(pnn, params);

  // prune  
  traln.clipNodeDefault(pn, pnn); 
  p->next->next->back = p->next->back = NULL; 

  path.append(bisectedBranch);
  path.append(BranchPlain(pn->number, pnn->number)); 
  nodeptr currentNode = rand.drawRandDouble01() ? pn : pnn;
  bool accepted = false; 
  while(not accepted )
    {
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      path.pushToStackIfNovel(BranchPlain(currentNode->number, n->number),traln); 
      currentNode = n;       
      accepted = rand.drawRandDouble01() < stopProb && path.size() > 2 ; 
    }

  // reset
  traln.clipNode(p->next, pn); 
  traln.setBranch(pnBranch, params);
  traln.clipNode(p->next->next, pnn); 
  traln.setBranch(pnnBranch, params);

  // a correction is necessary 
  if(path.at(2).hasNode(path.at(1).getPrimNode() ))
    path.at(1).setSecNode(p->number); 
  else if (path.at(2).hasNode(path.at(1).getSecNode() ))
    path.at(1).setPrimNode(p->number); 
  else 
    assert(0); 

  // for reasons of resetting the first branch in the path must be
  // identifyable later
  path.at(0) = bisectedBranch.getThirdBranch(traln, path.at(1)); 
}


BranchPlain ExtendedTBR::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
{
  auto canMove = [&](const BranchPlain &b) -> bool 
    { 
      auto result = b.isTipBranch(traln); 
      if( not result)
  	{
  	  auto desc = traln.getDescendents(b) ; 
  	  result = not(desc.first.isTipBranch(traln) && desc.second.isTipBranch(traln) ); 
  	  // tout << b << " with children " << std::get<0>(desc) << "," << std::get<1>(desc)<< std::endl ; 
  	}
      return result; 
    }; 

  auto bisectedBranch = BranchPlain{}; 
  auto movableA = bool{false}; 
  auto movableB = bool{false}; 

  while(not ( movableA || movableB ) )
    {
      bisectedBranch = TreeRandomizer::drawBranchWithInnerNode(traln,rand); 
      movableA = canMove(bisectedBranch); 
      movableB = canMove(bisectedBranch.getInverted()); 
    }

  return bisectedBranch; 
}


void ExtendedTBR::drawPaths(TreeAln &traln, Randomness &rand)
{
  auto modifiedPath1 = Path{}; 
  auto modifiedPath2 = Path{}; 

  auto bisectedBranch = determinePrimeBranch(traln, rand); 

  // determine, if a true TBR move can be executed
  auto canMove = [&](const BranchPlain &b) -> bool 
    { 
      if(traln.isTipNode(b.getPrimNode()))
	return false; 
      else 
	{
	  auto desc = traln.getDescendents(b) ; 
	  return not(desc.first.isTipBranch(traln) && desc.second.isTipBranch(traln)) ; 
	}
    }; 
  auto oneMovable = canMove(bisectedBranch);
  auto otherMovable = canMove(bisectedBranch.getInverted()); 

#ifdef DEBUG_TBR
  cout << "bisected branch is " << bisectedBranch << endl; 
#endif
  
  auto descOne = BranchPlain{}; 
  if(oneMovable)
    {
      buildPath(modifiedPath1, bisectedBranch, traln, rand); 
      descOne = modifiedPath1.at(modifiedPath1.size() -1 ); 
    }

  auto descOther = BranchPlain{}; 
  if(otherMovable)
    {
      
      buildPath(modifiedPath2, bisectedBranch.getInverted(), traln, rand); 
      descOther = modifiedPath2.at(modifiedPath2.size()-1); 
    }
  
  auto moveDescription = std::make_tuple( bisectedBranch, descOne, descOther); 
  _move.extractMoveInfo(traln, moveDescription , getSecondaryParameterView()); 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorBelief& prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{ 
  drawPaths(traln,rand);
  // TODO replace by absence of prior 

  // tout << "move is " << move << std::endl; 
  
  bool modifiesBl = false;   
  for(auto &v : _secondaryParameters)
    modifiesBl |= v->getCategory() == Category::BRANCH_LENGTHS; 

#ifdef NO_SEC_BL_MULTI
  modifiesBl = false; 
#endif

  if(modifiesBl)
    {
      // auto brPr = secondaryParameters.at(0)->getPrior();
      // move.multiplyBranches(traln, rand, hastings, prior,  multiplier, {brPr}); 
      assert(0); 
    }

  _move.applyToTree(traln, getSecondaryParameterView() );
}


void ExtendedTBR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln& traln, const BranchPlain &branchSuggestion)
{ 
  auto toEval = _move.getEvalBranch(traln);

  assert(toEval.exists(traln)); 
  
#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL " << toEval << std::endl; 
#endif

  auto dirtyNodes = _move.getDirtyNodes(traln, false);
  
  // tout << "dirtyNodes: " << dirtyNodes << std::endl; 

  for(auto &elem : dirtyNodes)
    evaluator.markDirty(traln,elem); 

  evaluator.evaluate(traln,toEval,false); 
}


void ExtendedTBR::resetState(TreeAln &traln)
{
  _move.revertTree(traln, getSecondaryParameterView() );
}



AbstractProposal* ExtendedTBR::clone() const 
{
  return new ExtendedTBR( *this );
}


std::vector<nat> ExtendedTBR::getInvalidatedNodes(const TreeAln& traln) const
{
  return _move.getDirtyNodes(traln, false);
} 
