#include "ExtendedTBR.hpp"
#include "output.h"
#include "branch.h"
#include "eval.h"

// #define EXPENSIVE_ETBR_VERIFY
// #define DEBUG_TBR


ExtendedTBR::ExtendedTBR( double _extensionProb, double _multiplier)
  :  extensionProbability(_extensionProb)
  , multiplier(_multiplier)
{
  name = "eTBR"; 
  category = TOPOLOGY; 
  relativeWeight = 5.;
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

  // prune  
  traln.clipNodeDefault(pn, pnn); 
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
  traln.clipNode(p->next, pn, z1); 
  traln.clipNode(p->next->next, pnn, z2); 

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

  branch bisectedBranch = getThirdBranch(tr, path.at(0), path.at(1)); 
  nodeptr disconnectedPtr = findNodeFromBranch(tr, bisectedBranch); 
  // cout << "bisected is " << bisectedBranch << "empty ptr is " << disconnectedPtr->number << endl;  

  branch b1 = constructBranch(path.getNthNodeInPath(0), path.getNthNodeInPath(1)); 
  // cout << "trying to find " << b1  <<endl; 
  nodeptr p1 = findNodeFromBranch(tr, b1); 

  int newConnect = path.getNthNodeInPath(2) ; 
  nodeptr pConn1 = findNodeFromBranch( tr, constructBranch(newConnect, getOtherNode(newConnect, path.at(1)))) ;
  double savedZ = pConn1->z[0]; 
  traln.clipNode(pConn1, p1, p1->z[0]);
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

  traln.clipNode(disconnectedPtr->next, iPtr, iPtr->z[0]); 
  traln.clipNode(disconnectedPtr->next->next, iPtr2, savedZ); 
}




void ExtendedTBR::executeTBR(TreeAln & traln)
{  
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 

#ifdef EXPENSIVE_ETBR_VERIFY
  double tlBefore = traln.getTreeLength(); 
#endif

  // prune the bisected branch from both sides, leaving behind two
  // unconnected trees
  pruneAndInsertOneSide(modifiedPath1, traln); 
  pruneAndInsertOneSide(modifiedPath2, traln); 

#ifdef EXPENSIVE_ETBR_VERIFY
  double tlAfter = traln.getTreeLength();
  if(fabs (tlAfter - tlBefore) > 1e-6   )
    {
      cout << "tbr changed bls where it should not" << endl; 
      assert(0); 
    }
#endif  
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorBelief& prior, double &hastings, Randomness &rand)
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
  assert(traln.getNumBranches() == 1); 
  // TODO replace by absence of prior 
  
  bool modifiesBl = false;   
  for(auto v : secVar)
    modifiesBl |= v.getCategory() == BRANCH_LENGTHS; 

  if(modifiesBl)
    {
      // find the branch length prior : very nasty 
      
      auto brPr = secVar.at(0).getPrior();
      for(int i = 0 ;i < modifiedPath1.size(); ++i)
	{
	  branch &b = modifiedPath1.at(i);  
	  modifiedPath1.multiplyBranch(traln, rand, b, multiplier, hastings, prior, brPr);       
	}
      for(int i = 0; i < modifiedPath2.size(); ++i)
	{
	  branch &b = modifiedPath2.at(i); 
	  modifiedPath2.multiplyBranch(traln, rand, b, multiplier, hastings, prior, brPr); 
	}

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



static void disorientHelper(TreeAln &traln, nodeptr p)
{
  if(NOT traln.isTipNode(p))
    {
      p->next->x = p->next->next->x = 0; 
      p->x = 1 ; 
      // cout << "disorienting " << p->number << endl; 
    }
  // else 
    // cout << "NOT disorienting " << p->number << endl; 
}


static void destroyOrientationAlongPath(TreeAln &traln, Path &path)
{
  tree *tr = traln.getTr(); 
  // TODO necessary for the last node? 
  for(int i = 0 ; i < path.getNumberOfNodes()-1; ++i)
    {
      nodeptr p = findNodeFromBranch(tr,  constructBranch(path.getNthNodeInPath(i) ,path.getNthNodeInPath(i+1) ));
      disorientHelper(traln, p); 
    }
  
  // nodeptr p = findNodeFromBranch(tr , constructBranch(path.getNthNodeInPath(path.getNumberOfNodes()-1), 
  // 						      path.getNthNodeInPath(path.getNumberOfNodes()-2))); 
  // disorientHelper(traln, p); 
}



static void evalHelper(Path &path, TreeAln& traln )
{  
  tree *tr = traln.getTr();
  nodeptr p = findNodeFromBranch ( tr, constructBranch(path.getNthNodeInPath(path.getNumberOfNodes()-1), 
						       path.getNthNodeInPath(1))); 		       
  newViewGenericWrapper(traln, p, FALSE);  

  Path pathCopy(path); 
  pathCopy.pop();   
  branch prunedBranch = constructBranch(pathCopy.getNthNodeInPath(0), pathCopy.getNthNodeInPath(2)); 
  pathCopy.popFront();
  pathCopy.at(0) = prunedBranch; 
  pathCopy.append(constructBranch(path.getNthNodeInPath(1),path.getNthNodeInPath(path.getNumberOfNodes()-2))); 
  pathCopy.reverse(); 
  // cout << "path for eval is " << pathCopy << endl; 
  destroyOrientationAlongPath(traln, pathCopy); 
  p = findNodeFromBranch(tr, constructBranch( pathCopy.getNthNodeInPath(0), pathCopy.getNthNodeInPath(1))); 
  if(NOT traln.isTipNode(p))
    newViewGenericWrapper(traln, p,FALSE); 
}



void ExtendedTBR::evaluateProposal(TreeAln& traln, PriorBelief& prior)
{  
#if 1 
  tree *tr = traln.getTr(); 

  branch bisectedBranch = constructBranch(modifiedPath1.getNthNodeInPath(1), modifiedPath2.getNthNodeInPath(1)); 
  
  // evaluate both ends of the paths 
  evalHelper(modifiedPath1, traln); 
  evalHelper(modifiedPath2, traln); 
  nodeptr p = findNodeFromBranch(tr, bisectedBranch);
  evaluateGenericWrapper(traln, p, FALSE);
  
#else 

  // IMPORTANT TODO !!! 
#ifdef DEBUG_TBR
  cout << "lnl before = "  << traln.getTr()->likelihood << endl; 
#endif
  evaluateGenericWrapper(traln,traln.getTr()->start, TRUE );
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
  
  int lastNode = path.getNthNodeInPath(path.getNumberOfNodes() -1 ),
    s2lastNode =  path.getNthNodeInPath(path.getNumberOfNodes()-2);   
  int subTreeNode = path.at(0).thisNode; 

  nodeptr p1 = findNodeFromBranch(tr, constructBranch(lastNode, subTreeNode)),
    p2 = findNodeFromBranch(tr, constructBranch(s2lastNode, subTreeNode));     

  // double defaultVal = 0.123; 

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
  traln.clipNodeDefault(p1,p2); 
  
  // reinsert the subtree 
  branch insertBranch = constructBranch(path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 
  nodeptr iPtr = findNodeFromBranch(tr, insertBranch),
    iPtr2 = findNodeFromBranch(tr, invertBranch(insertBranch)); 
  
  traln.clipNodeDefault(iPtr, subTreePtr->next); 
  traln.clipNodeDefault(iPtr2, subTreePtr->next->next); 
}


void ExtendedTBR::resetState(TreeAln &traln, PriorBelief& prior)
{
  resetOneSide(traln, modifiedPath1); 
  resetOneSide(traln, modifiedPath2); 

  modifiedPath1.restoreBranchLengthsPath(traln, prior); 
  modifiedPath2.restoreBranchLengthsPath(traln, prior); 


#ifdef EFFICIENT
  assert(0); 
#endif
 
  debug_checkTreeConsistency(traln.getTr()); 
  debug_printTree(traln);   
}



AbstractProposal* ExtendedTBR::clone() const 
{
  return new ExtendedTBR( *this );
// extensionProbability, multiplier
}
