#include "ExtendedTBR.hpp"
#include "output.h"
#include "eval.h"


ExtendedTBR::ExtendedTBR( double _extensionProb, double _multiplier)
  :  extensionProbability(_extensionProb)
  , multiplier(_multiplier)
{
  name = "eTBR"; 
  category = TOPOLOGY; 
  relativeWeight = 5.;
}


static void buildPath(Path &path, Branch bisectedBranch, TreeAln &traln, Randomness &rand, double stopProb )
{
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
  tree *tr  = traln.getTr();

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
  while(NOT accepted )
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
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ) ; 
  tree *tr = traln.getTr();

  Branch bisectedBranch = {0,0}; 
  

  assert(tr->mxtips > 8 ); 

  nodeptr
    p1, p2;  

  do
    {
      Branch bisectedBranch  = traln.drawInnerBranchUniform(rand); 

      p1 = bisectedBranch.findNodePtr(traln ); 
      p2 = bisectedBranch.getInverted().findNodePtr(traln); 

      assert(NOT isTip(p1->number, tr->mxtips) && NOT isTip(p2->number, tr->mxtips)); 

    } while( (traln.isTipNode(p1->next->back) && traln.isTipNode(p1->next->next->back) ) 
	     || (traln.isTipNode(p2->next->back)  && traln.isTipNode(p2->next->next->back) ) ) ; 

				// variables otherwise

#ifdef DEBUG_TBR
  cout << "bisected branch is " << bisectedBranch << endl; 
#endif

  buildPath(modifiedPath1, bisectedBranch, traln, rand, extensionProbability); 
  buildPath(modifiedPath2, bisectedBranch.getInverted(), traln, rand, extensionProbability); 
}



static void pruneAndInsertOneSide(Path &path, TreeAln &traln)
{   
  tree *tr  = traln.getTr(); 

  Branch bisectedBranch = path.at(0).getThirdBranch(traln , path.at(1)); 
  nodeptr disconnectedPtr = bisectedBranch.findNodePtr(traln ); 

  Branch b1(path.getNthNodeInPath(0), path.getNthNodeInPath(1)); 
  nodeptr p1 = b1.findNodePtr(traln ); 

  int newConnect = path.getNthNodeInPath(2) ; 
  nodeptr pConn1 = Branch(newConnect, path.at(1).getOtherNode(newConnect )).findNodePtr( traln ) ;
  double savedZ = pConn1->z[0]; 
  traln.clipNode(pConn1, p1, p1->z[0]);
  disconnectedPtr->next->back =  disconnectedPtr->next->next->back  = NULL; 
#ifdef DEBUG_TBR
  cout << "hooking " << pConn1->number << "," <<p1->number << "(" << p1->z[0] << ")" << endl; 
#endif

  Branch insertBranch(path.getNthNodeInPath(path.getNumberOfNodes( )- 1 ), path.getNthNodeInPath(path.getNumberOfNodes()-2) );   
#ifdef DEBUG_TBR
  cout << "insert branch is " << insertBranch << endl; 
#endif
  nodeptr iPtr = insertBranch.findNodePtr(traln ),
    iPtr2 = insertBranch.getInverted().findNodePtr(traln); 

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

  assert(traln.getNumBranches() == 1); 
  // TODO replace by absence of prior 
  
  bool modifiesBl = false;   
  for(auto v : secVar)
    modifiesBl |= v.getCategory() == BRANCH_LENGTHS; 

#ifdef NO_SEC_BL_MULTI
  modifiesBl = false; 
#endif

  if(modifiesBl)
    {
      // find the branch length prior : very nasty 
      
      auto brPr = secVar.at(0).getPrior();
      for(nat i = 0 ;i < modifiedPath1.size(); ++i)
	{
	  Branch &b = modifiedPath1.at(i);  
	  modifiedPath1.multiplyBranch(traln, rand, b, multiplier, hastings, prior, brPr);       
	}
      for(nat i = 0; i < modifiedPath2.size(); ++i)
	{
	  Branch &b = modifiedPath2.at(i); 
	  modifiedPath2.multiplyBranch(traln, rand, b, multiplier, hastings, prior, brPr); 
	}

    }

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
    }
}


static void destroyOrientationAlongPath(TreeAln &traln, Path &path)
{
  tree *tr = traln.getTr(); 
  // TODO necessary for the last node? 
  for(int i = 0 ; i < path.getNumberOfNodes()-1; ++i)
    {
      nodeptr p = Branch(path.getNthNodeInPath(i) ,path.getNthNodeInPath(i+1) ).findNodePtr(traln );
      disorientHelper(traln, p); 
    }
  
  // nodeptr p = findNodeFromBranch(tr , constructBranch(path.getNthNodeInPath(path.getNumberOfNodes()-1), 
  // 						      path.getNthNodeInPath(path.getNumberOfNodes()-2))); 
  // disorientHelper(traln, p); 
}



static void evalHelper(Path &path, TreeAln& traln )
{  
  tree *tr = traln.getTr();
  nodeptr p = Branch(path.getNthNodeInPath(path.getNumberOfNodes()-1), path.getNthNodeInPath(1)).findNodePtr ( traln ); 		       
  newViewGenericWrapper(traln, p, FALSE);  

  Path pathCopy(path); 
  pathCopy.pop();   
  Branch prunedBranch(pathCopy.getNthNodeInPath(0), pathCopy.getNthNodeInPath(2)); 
  pathCopy.popFront();
  pathCopy.at(0) = prunedBranch; 
  pathCopy.append(Branch(path.getNthNodeInPath(1),path.getNthNodeInPath(path.getNumberOfNodes()-2))); 
  pathCopy.reverse(); 
  // cout << "path for eval is " << pathCopy << endl; 
  destroyOrientationAlongPath(traln, pathCopy); 
  p = Branch( pathCopy.getNthNodeInPath(0), pathCopy.getNthNodeInPath(1)).findNodePtr(traln ); 
  if(NOT traln.isTipNode(p))
    newViewGenericWrapper(traln, p,FALSE); 
}



void ExtendedTBR::evaluateProposal(TreeAln& traln, PriorBelief& prior)
{  
#if 1 
  tree *tr = traln.getTr(); 

  Branch bisectedBranch(modifiedPath1.getNthNodeInPath(1), modifiedPath2.getNthNodeInPath(1)); 
  
  // evaluate both ends of the paths 
  evalHelper(modifiedPath1, traln); 
  evalHelper(modifiedPath2, traln); 
  nodeptr p = bisectedBranch.findNodePtr(traln );
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
  int subTreeNode = path.at(0).getPrimNode(); 

  nodeptr p1 = Branch(lastNode, subTreeNode).findNodePtr(traln ),
    p2 = Branch(s2lastNode, subTreeNode).findNodePtr(traln );     

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
  Branch insertBranch(path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 
  nodeptr iPtr = insertBranch.findNodePtr(traln ),
    iPtr2 = insertBranch.getInverted().findNodePtr(traln); 
  
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
 
  debug_checkTreeConsistency(traln); 
  debug_printTree(traln);   
}



AbstractProposal* ExtendedTBR::clone() const 
{
  return new ExtendedTBR( *this );
// extensionProbability, multiplier
}
