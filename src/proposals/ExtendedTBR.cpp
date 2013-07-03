#include "ExtendedTBR.hpp"


// TODO the disorient is still  very inefficient 

ExtendedTBR::ExtendedTBR( double _extensionProb, double _multiplier)
  :  extensionProbability(_extensionProb)
  , multiplier(_multiplier)
{
  name = "eTBR"; 
  category = Category::TOPOLOGY; 
  relativeWeight = 5.;
}


void ExtendedTBR::buildPath(Path &path, Branch bisectedBranch, TreeAln &traln, Randomness &rand )
{
  double stopProb = extensionProbability; 

  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 

  nodeptr p= bisectedBranch.findNodePtr(traln ); 
  path.clear();

  nodeptr pn = p->next->back,
    pnn = p->next->next->back; 
  
  double z1 = pn->z[0],
    z2 = pnn->z[0]; 

  // prune  
  traln.clipNodeDefault(pn, pnn); 
  p->next->next->back = p->next->back = NULL; 
  path.append(bisectedBranch);
  path.append(Branch(pn->number, pnn->number)); 
  nodeptr currentNode = rand.drawRandDouble01() ? pn : pnn;
  bool accepted = false; 
  while(not accepted )
    {
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      path.pushToStackIfNovel(Branch(currentNode->number, n->number),traln); 
      currentNode = n;       
      accepted = rand.drawRandDouble01() < stopProb && path.size() > 2 ; 
    }

  // reset
  traln.clipNode(p->next, pn, z1); 
  traln.clipNode(p->next->next, pnn, z2); 

  // a correction is necessary 
  if(path.at(2).nodeIsInBranch(path.at(1).getPrimNode() ))
    path.at(1).setSecNode(p->number); 
  else if (path.at(2).nodeIsInBranch(path.at(1).getSecNode() ))
    path.at(1).setPrimNode(p->number); 
  else 
    assert(0); 

  // for reasons of resetting the first branch in the path must be
  // identifyable later
  path.at(0) = bisectedBranch.getThirdBranch(traln, path.at(1)); 
}


void ExtendedTBR::drawPaths(TreeAln &traln, Randomness &rand)
{
  Path modifiedPath1, 
    modifiedPath2; 

  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ) ; 
  tree *tr = traln.getTr();

  Branch bisectedBranch = {0,0}; 
  

  assert(tr->mxtips > 8 ); 

  nodeptr
    p1, p2;  

  do
    {
      bisectedBranch  = traln.drawInnerBranchUniform(rand ); 

      p1 = bisectedBranch.findNodePtr(traln ); 
      p2 = bisectedBranch.getInverted().findNodePtr(traln); 

      assert(not isTip(p1->number, tr->mxtips) && not isTip(p2->number, tr->mxtips)); 

    } while( (traln.isTipNode(p1->next->back) && traln.isTipNode(p1->next->next->back) ) 
	     || (traln.isTipNode(p2->next->back)  && traln.isTipNode(p2->next->next->back) ) ) ; 

  // variables otherwise

#ifdef DEBUG_TBR
  cout << "bisected branch is " << bisectedBranch << endl; 
#endif

  buildPath(modifiedPath1, bisectedBranch, traln, rand); 
  buildPath(modifiedPath2, bisectedBranch.getInverted(), traln, rand); 

  move.extractMoveInfo(traln, {			 
      bisectedBranch, modifiedPath1.at(modifiedPath1.size() -1 ),
	bisectedBranch.getInverted(), modifiedPath2.at(modifiedPath2.size()-1) } ); 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorBelief& prior, double &hastings, Randomness &rand)
{ 
  drawPaths(traln,rand);

  assert(traln.getNumBranches() == 1); 
  // TODO replace by absence of prior 
  
  bool modifiesBl = false;   
  for(auto &v : secVar)
    modifiesBl |= v->getCategory() == Category::BRANCH_LENGTHS; 

#ifdef NO_SEC_BL_MULTI
  modifiesBl = false; 
#endif

  if(modifiesBl)
    {
      auto brPr = secVar.at(0)->getPrior();
      move.multiplyBranches(traln, rand, hastings, prior,  multiplier, {brPr}); 
    }

  move.applyToTree(traln);
}


void ExtendedTBR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln& traln, PriorBelief& prior)
{ 
  Branch toEval = move.getEvalBranch(traln);
  auto p = toEval.findNodePtr(traln); 
  move.disorientAtNode(traln,p->back);     

  evaluator.evaluate(traln,toEval,false); 
}


void ExtendedTBR::resetState(TreeAln &traln, PriorBelief& prior)
{
  move.revertTree(traln,prior);

#ifdef EFFICIENT
  assert(0); 
#endif

  // debug_checkTreeConsistency(traln); 
  // debug_printTree(traln);   
}



AbstractProposal* ExtendedTBR::clone() const 
{
  return new ExtendedTBR( *this );
}

