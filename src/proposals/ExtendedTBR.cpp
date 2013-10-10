#include "ExtendedTBR.hpp"
#include "Category.hpp"
#include "TreeRandomizer.hpp"


// TODO the disorient is still  very inefficient 

ExtendedTBR::ExtendedTBR( double _extensionProb, double _multiplier)
  : AbstractProposal(Category::TOPOLOGY , "eTBR")
  , extensionProbability(_extensionProb)
  , multiplier(_multiplier)
{
  relativeWeight = 5.;
  needsFullTraversal = false;
}


void ExtendedTBR::buildPath(Path &path, BranchPlain bisectedBranch, TreeAln &traln, Randomness &rand )
{
  auto params = getBranchLengthsParameterView();

  double stopProb = extensionProbability; 

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

  nodeptr
    p1, p2;  

  auto bisectedBranch = BranchPlain(); 
  do
    {
      bisectedBranch = TreeRandomizer::drawInnerBranchUniform(traln, rand ); 

      p1 = bisectedBranch.findNodePtr(traln ); 
      p2 = bisectedBranch.getInverted().findNodePtr(traln); 

      assert(not traln.isTipNode(p1) && not traln.isTipNode(p2)); 

    } while( (traln.isTipNode(p1->next->back) && traln.isTipNode(p1->next->next->back) ) 
	     || (traln.isTipNode(p2->next->back)  && traln.isTipNode(p2->next->next->back) ) ) ;   
  
  return bisectedBranch; 
}


void ExtendedTBR::drawPaths(TreeAln &traln, Randomness &rand)
{
  Path modifiedPath1, 
    modifiedPath2; 

  // int numBranches = traln.getNumBranches(); 
  // assert(numBranches == 1 ) ; 
  tree *tr = traln.getTr();

  auto bisectedBranch = BranchPlain{0,0}; 
  assert(tr->mxtips > 8 ); 
  
  bisectedBranch = determinePrimeBranch(traln, rand); 
  // auto p1 = bisectedBranch.findNodePtr(traln); 
  // auto p2 = bisectedBranch.getInverted().findNodePtr(traln); 

  // variables otherwise

#ifdef DEBUG_TBR
  cout << "bisected branch is " << bisectedBranch << endl; 
#endif

  buildPath(modifiedPath1, bisectedBranch, traln, rand); 
  buildPath(modifiedPath2, bisectedBranch.getInverted(), traln, rand); 

  move.extractMoveInfo(traln, {			 
      bisectedBranch, modifiedPath1.at(modifiedPath1.size() -1 ),
	bisectedBranch.getInverted(), modifiedPath2.at(modifiedPath2.size()-1) } , getSecondaryParameterView()); 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorBelief& prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{ 
  drawPaths(traln,rand);

  // TODO replace by absence of prior 
  
  bool modifiesBl = false;   
  for(auto &v : secondaryParameters)
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

  move.applyToTree(traln, getSecondaryParameterView() );
}


void ExtendedTBR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln& traln, const BranchPlain &branchSuggestion)
{ 
  auto toEval = move.getEvalBranch(traln);
  
#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL " << toEval << std::endl; 
#endif

  auto dirtyNodes = move.getDirtyNodes();
  
  for(auto &elem : dirtyNodes)
    evaluator.markDirty(traln,elem); 

  evaluator.evaluate(traln,toEval,false); 
}


void ExtendedTBR::resetState(TreeAln &traln)
{
  move.revertTree(traln, getSecondaryParameterView() );
}



AbstractProposal* ExtendedTBR::clone() const 
{
  return new ExtendedTBR( *this );
}


std::vector<nat> ExtendedTBR::getInvalidatedNodes(const TreeAln& traln) const
{
  return move.getDirtyNodes();
} 
