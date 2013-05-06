#include "ExtendedTBR.hpp"
#include "output.h"
#include "branch.h"
#include "topology-utils.h"


// #define DEBUG_TBR


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


static void buildPath(Path &path, branch bisectedBranch, TreeAln &traln, Randomness &rand, double stopProb )
{
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
  tree *tr  = traln.getTr();

  nodeptr p=  findNodeFromBranch(tr, bisectedBranch); 
  path.clear();

  nodeptr pn = p->next->back,
    pnn = p->next->next->back; 
  
  double z1 = pn->z[0],
    z2 = pnn->z[0]; 

  double someDefault = 0.123; 

  // prune  
  hookup(pn, pnn, &someDefault, numBranches); 
  p->next->next->back = p->next->back = NULL; 
  path.append(bisectedBranch);
  path.append(constructBranch(pn->number, pnn->number)); 
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

  // for reasons of resetting the first branch in the path must be
  // identifyable later
  path.at(0) = getThirdBranch(tr, bisectedBranch, path.at(1)); 
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

#ifdef DEBUG_TBR
  cout << "bisected branch is " << bisectedBranch << endl; 
#endif

  buildPath(modifiedPath1, bisectedBranch, traln, rand, extensionProbability); 
  buildPath(modifiedPath2, invertBranch(bisectedBranch), traln, rand, extensionProbability); 
}



static void pruneAndInsertOneSide(Path &path, TreeAln &traln)
{   
  tree *tr  = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 

  branch bisectedBranch = getThirdBranch(tr, path.at(0), path.at(1)); 
  nodeptr disconnectedPtr = findNodeFromBranch(tr, bisectedBranch); 
  // cout << "bisected is " << bisectedBranch << "empty ptr is " << disconnectedPtr->number << endl;  

  branch b1 = constructBranch(path.getNthNodeInPath(0), path.getNthNodeInPath(1)); 
  // cout << "trying to find " << b1  <<endl; 
  nodeptr p1 = findNodeFromBranch(tr, b1); 

  int newConnect = path.getNthNodeInPath(2) ; 
  nodeptr pConn1 = findNodeFromBranch( tr, constructBranch(newConnect, getOtherNode(newConnect, path.at(1)))) ;
  double savedZ = pConn1->z[0]; 
  hookup(pConn1, p1, p1->z, numBranches);
  disconnectedPtr->next->back =  disconnectedPtr->next->next->back  = NULL; 
#ifdef DEBUG_TBR
  cout << "hooking " << pConn1->number << "," <<p1->number << "(" << p1->z[0] << ")" << endl; 
#endif

  branch insertBranch  = constructBranch(path.getNthNodeInPath(path.getNumberOfNodes( )- 1 ), path.getNthNodeInPath(path.getNumberOfNodes()-2) );   
#ifdef DEBUG_TBR
  cout << "insert branch is " << insertBranch << endl; 
#endif
  nodeptr iPtr = findNodeFromBranch(tr, insertBranch),
    iPtr2 = findNodeFromBranch(tr, invertBranch(insertBranch)); 


  hookup(disconnectedPtr->next, iPtr, iPtr->z, numBranches); 
  hookup(disconnectedPtr->next->next, iPtr2, &savedZ, numBranches); 
}


void ExtendedTBR::executeTBR(TreeAln & traln)
{  
  // tree *tr = traln.getTr();
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 

  // prune the bisected branch from both sides, leaving behind two
  // unconnected trees
  pruneAndInsertOneSide(modifiedPath1, traln); 
  pruneAndInsertOneSide(modifiedPath2, traln); 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorManager& prior, double &hastings, Randomness &rand)
{
  // tree *tr = traln.getTr(); 
#ifdef DEBUG_TBR
  cout << "before: " << traln << endl; 
#endif

  modifiedPath1.clear();
  modifiedPath2.clear();

  drawPaths(traln,rand);

  modifiedPath1.saveBranchLengthsPath(traln); 
  modifiedPath2.saveBranchLengthsPath(traln); 

#ifdef TBR_MULTIPLY_BL  
  for(int i = 0 ;i < modifiedPath1.size(); ++i)
    {
      branch &b = modifiedPath1.at(i);  
      modifiedPath1.multiplyBranch(traln, rand, b, multiplier, hastings);       
    }
  for(int i = 0; i < modifiedPath2.size(); ++i)
    {
      branch &b = modifiedPath2.at(i); 
      modifiedPath2.multiplyBranch(traln, rand, b, multiplier, hastings); 
    }
#endif
  
#ifdef DEBUG_TBR
  cout << modifiedPath1 <<  endl; 
  cout << modifiedPath2 <<  endl; 
#endif

  executeTBR(traln);
  
#ifdef DEBUG_TBR
  cout << "after: " << traln << endl; 
#endif
}


void ExtendedTBR::evaluateProposal(TreeAln& traln, PriorManager& prior)
{
#if 1 
  assert(0); 



#else 

  // IMPORTANT TODO !!! 
#ifdef DEBUG_TBR
  cout << "lnl before = "  << traln.getTr()->likelihood << endl; 
#endif
  evaluateGenericWrapper(chain,traln.getTr()->start, TRUE );
#ifdef DEBUG_TBR
  cout << "lnl after = " << traln.getTr()->likelihood << endl; 
#endif
#endif
}




/** 
    @notice dont care so much about branch lengths here, we have a
    reset later anyway.
 */ 
static void resetOneSide(TreeAln &traln, Path& path)
{  
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 
  
  int lastNode = path.getNthNodeInPath(path.getNumberOfNodes() -1 ),
    s2lastNode =  path.getNthNodeInPath(path.getNumberOfNodes()-2);   
  int subTreeNode = path.at(0).thisNode; 

  nodeptr p1 = findNodeFromBranch(tr, constructBranch(lastNode, subTreeNode)),
    p2 = findNodeFromBranch(tr, constructBranch(s2lastNode, subTreeNode));     

  double defaultVal = 0.123; 

  // correctly orient the subtree ptr 
  nodeptr subTreePtr = NULL; 
  {
    subTreePtr = p1->back; 
    p1->back->back = NULL; 
    p2->back->back = NULL;
    while(subTreePtr->back == NULL)
      subTreePtr = subTreePtr->next;     
  }
  
  // prune the subtree  
  assert(p1->back->number == p2->back->number); 
  hookup(p1,p2,&defaultVal, numBranches); 
  
  // reinsert the subtree 
  branch insertBranch = constructBranch(path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 
  nodeptr iPtr = findNodeFromBranch(tr, insertBranch),
    iPtr2 = findNodeFromBranch(tr, invertBranch(insertBranch)); 
  
  hookup(iPtr, subTreePtr->next, &defaultVal,numBranches); 
  hookup(iPtr2, subTreePtr->next->next, &defaultVal, numBranches); 
}


void ExtendedTBR::resetState(TreeAln &traln, PriorManager& prior)
{
  resetOneSide(traln, modifiedPath1); 
  resetOneSide(traln, modifiedPath2); 

  modifiedPath1.restoreBranchLengthsPath(traln); 
  modifiedPath2.restoreBranchLengthsPath(traln); 

  debug_checkTreeConsistency(traln.getTr()); 
  debug_printTree(traln);   
}
