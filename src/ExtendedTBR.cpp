#include "ExtendedTBR.hpp"
#include "branch.h"




static void buildPath(Path &path, branch bisectedBranch, TreeAln &traln, Randomness &rand, double stopProb )
{
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
  tree *tr  = traln.getTr();

  nodeptr p=  findNodeFromBranch(tr, bisectedBranch); 
  path.clearStack();

  nodeptr pn = p->next->back,
    pnn = p->next->next->back; 
  
  double z1 = pn->z[0],
    z2 = pnn->z[0]; 

  double someDefault = 0.123; 

  // prune  
  hookup(pn, pnn, &someDefault, numBranches); 
  p->next->next->back = p->next->back = NULL; 
  path.pushStack(bisectedBranch);
  path.pushStack(constructBranch(pn->number, pnn->number)); 
  nodeptr currentNode = rand.drawRandDouble01() ? pn : pnn;
  bool accepted = false; 
  while(NOT accepted )
    {
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      path.pushToStackIfNovel(constructBranch(currentNode->number, n->number),tr->mxtips); 
      currentNode = n;       
      accepted = rand.drawRandDouble01() < stopProb && path.size() > 2 ; 
    }

  // reset
  hookup(p->next, pn, &z1, numBranches); 
  hookup(p->next->next, pnn, &z2, numBranches); 

  // a correction is necessary 
  if(nodeIsInBranch(path.at(1).thisNode, path.at(2)))
    path.at(1).thatNode = p->number; 
  else if (nodeIsInBranch(path.at(1).thatNode, path.at(2)))
    path.at(1).thisNode = p->number;     
  else 
    assert(0); 
}



ExtendedTBR::ExtendedTBR( Chain *_chain, double _relativeProbability, double _extensionProb, double _multiplier)
  : chain(_chain)
  ,  extensionProbability(_extensionProb)
  , multiplier(_multiplier)
{
  name = "eTBR"; 
  relativeProbability = _relativeProbability; 
  category = TOPOLOGY; 
  ptype = E_TBR; 
}



void ExtendedTBR::autotune()
{
  // cannot tune that 
}


void ExtendedTBR::drawPaths(TreeAln &traln, Randomness &rand)
{
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ) ; 
  tree *tr = traln.getTr();

  branch bisectedBranch = {0,0}; 

  assert(tr->mxtips > 8 ); 

  nodeptr
    p1, p2;  

  do
    {
      bisectedBranch = rand.drawInnerBranchUniform(traln); 
      p1 = findNodeFromBranch(tr, bisectedBranch); 
      p2 = findNodeFromBranch(tr, invertBranch( bisectedBranch)); 

      assert(NOT isTip(p1->number, tr->mxtips) && NOT isTip(p2->number, tr->mxtips)); 
      
      // we do not want a subtree that is a tip-tip case 
      if(NOT (( isTip( p1->next->back->number, tr->mxtips) && isTip( p1->next->next->back->number, tr->mxtips) ) 
	      || ( isTip( p2->next->back->number, tr->mxtips) && isTip( p2->next->next->back->number, tr->mxtips) ) ))
	break; 
      
      
    } while( TRUE  ); 		// just because of the many nasty
				// variables otherwise

  
  buildPath(modifiedPath1, bisectedBranch, traln, rand, extensionProbability); 
  buildPath(modifiedPath2, invertBranch(bisectedBranch), traln, rand, extensionProbability); 
}



static void pruneOneSide(Path &path, TreeAln &traln, double &savedZ)
{
  tree *tr  = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 

  branch b1 = getThirdBranch(tr, path.at(0), path.at(1)); 
  nodeptr p1 = findNodeFromBranch(tr, b1)->back; 

  int newConnect = path.getNthNodeInPath(2) ; 
  nodeptr pConn1 = findNodeFromBranch( tr, constructBranch(newConnect, getOtherNode(newConnect, path.at(1)))) ;
  savedZ = pConn1->z[0]; 
  hookup(pConn1, p1, p1->z, numBranches);
  
  cout << "hooking " << pConn1->number << "," <<p1->number << "(" << p1->z[0] << ")" << endl; 
}



static void insertOneSide(Path &path, TreeAln &traln, double savedZ)
{  
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 

  int lastNode = path.getNthNodeInPath(path.getNumberOfNodes()-1),
    s2LastNode = path.getNthNodeInPath(path.getNumberOfNodes()-2); 

  nodeptr tmp = findNodeFromBranch(tr, path.at(0)),
    iPtr1 = tmp->next,
    iPtr2 = tmp->next->next; 

  nodeptr lastPtr = findNodeFromBranch(tr, constructBranch(lastNode, s2LastNode)),
    s2LastPtr = findNodeFromBranch(tr, constructBranch(s2LastNode, lastNode)); 

  hookup(iPtr2, lastPtr, lastPtr->z, numBranches); 
  hookup(iPtr1, s2LastPtr, &savedZ, numBranches); 

  cout << "hooked up " << iPtr1->number << "," << s2LastPtr->number << "(" << iPtr1->z[0] << ")\tand\t"  << iPtr2->number << "," << lastPtr->number << "(" << lastPtr->z[0] << ")" << endl; 
}



void ExtendedTBR::executeTBR(TreeAln & traln)
{  
  // tree *tr = traln.getTr();
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 

  assert(branchEqualUndirected(modifiedPath1.at(0), modifiedPath2.at(0))); 
  
  // prune the bisected branch from both sides, leaving behind two
  // unconnected trees
  double saved1 = 0,  saved2 = 0; 
  pruneOneSide(modifiedPath1, traln, saved1); 
  pruneOneSide(modifiedPath2, traln, saved2); 
  
  insertOneSide(modifiedPath1,traln, saved1); 
  insertOneSide(modifiedPath2, traln, saved2); 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorManager& prior, double &hastings, Randomness &rand)
{
  // tree *tr = traln.getTr(); 

  cout << "before: " << traln << endl; 

  modifiedPath1.clearStack();
  modifiedPath2.clearStack();

  drawPaths(traln,rand);

  modifiedPath1.saveBranchLengthsPath(traln); 
  modifiedPath2.saveBranchLengthsPath(traln); 

#ifdef TBR_MULTIPLY_BL  
  int ctr = 0; 
  for(auto  b : modifiedPath1.getStack() )
    {
      if(ctr == 0)continue; // do nothing to the first guy! -- while this could also be modified  	 
      modifiedPath1.multiplyBranch(traln, rand, b, multiplier, hastings);       
    }
  ctr = 0; 
  for(auto b : modifiedPath2.getStack())
    {
      if(ctr == 0) continue; 
      modifiedPath2.multiplyBranch(traln, rand, b, multiplier, hastings); 
    }
#endif

  cout << modifiedPath1 <<  endl; 
  cout << modifiedPath2 <<  endl; 

  executeTBR(traln);

  cout << "after: " << traln << endl; 

}


void ExtendedTBR::evaluateProposal(TreeAln& traln, PriorManager& prior)
{
  // IMPORTANT TODO !!! 
  cout << "lnl before = "  << traln.getTr()->likelihood << endl; 
  evaluateGenericWrapper(chain,traln.getTr()->start, TRUE );
  cout << "lnl after = " << traln.getTr()->likelihood << endl; 
}


void ExtendedTBR::resetState(TreeAln &traln, PriorManager& prior)
{
  
}

