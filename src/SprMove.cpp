 #include <set>

#include "SprMove.hpp"
#include "branch.h"

// TODO make branch stuff more optional (performance)


void SprMove::applyToTree(TreeAln &traln) const
{
  applyPathAsESPR(traln, path); 
}

void SprMove::revertTree(TreeAln &traln, PriorBelief &prior) const
{
  resetAlongPathForESPR(traln, prior, path) ;   
  path.restoreBranchLengthsPath(traln, prior); 
}

 
void SprMove::disorientAtNode(TreeAln &traln, nodeptr p) const
{
  int first = (int)path.getNthNodeInPath(0),
    last = (int)path.getNthNodeInPath(path.getNumberOfNodes()-1); 

  if(NOT path.nodeIsOnPath(p->number) || traln.isTipNode(p) || p->number == first || p->number == last )
    return; 

  disorientHelper(traln, p);
  
  disorientAtNode( traln, p->next->back); 
  disorientAtNode( traln, p->next->next->back);
}


void SprMove::extractMoveInfo(const TreeAln &traln, vector<Branch> description)
{
#ifdef EFFICIENT
  // looks very clunky 
  assert(0); 
#endif

  Branch chosen = description.at(1); 
  Branch prunedTree = description.at(0); 

  nodeptr p = prunedTree.findNodePtr(traln);
  
  path.clear();   
  path.findPath(traln, p, chosen.findNodePtr(traln));
  path.reverse();   

  branch tmp = path.at(0); 
  Branch Tmp = Branch(tmp.thisNode, tmp.thatNode); 

  nodeptr aNode = Tmp.findNodePtr(traln); 
  if(traln.isTipNode(aNode))
    aNode = findNodeFromBranch(traln.getTr(), invertBranch(path.at(0)));
  assert(not traln.isTipNode(aNode)); 
  branch b = path.at(0); 
  while(branchEqualUndirected(b, path.at(0)) || branchEqualUndirected(b,constructBranch(p->number, p->back->number))) 
    {
      aNode = aNode->next; 
      b = constructBranch(aNode->number, aNode->back->number) ;
    }
  path.reverse();
  path.append(b);  
  path.reverse();

  path.saveBranchLengthsPath(traln); 
  // =/ 
} 


AbstractMove* SprMove::clone() const
{
  return new SprMove; 
} 




Branch SprMove::getEvalBranch(const TreeAln &traln) const
{  
  // LEGACY 
  auto *tr = traln.getTr(); 
  
  branch b = getThirdBranch(tr, path.at(0), path.at(1)); 
  branch lastBranch = path.at(path.size()-1); 

  branch bA = constructBranch(b.thisNode, lastBranch.thisNode), 
    bB = constructBranch(b.thisNode,lastBranch.thatNode); 

  assert(branchExists(tr, bA) && branchExists(tr,bB)); 

  branch futureRoot = getThirdBranch(tr, bA, bB ); 
  
  return Branch(futureRoot.thisNode, futureRoot.thatNode); 
}



void SprMove::multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, shared_ptr<AbstractPrior> brPr)  const 
{  
  multiplyAlongBranchESPR(traln, rand, hastings, prior,  path, multiplier,  brPr); 
} 


/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void SprMove::multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior,  const Path &modifiedPath, double multiplier, shared_ptr<AbstractPrior> brPr) const
{
  assert(modifiedPath.size() >= 2); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
 
  int sTNode = modifiedPath.getNthNodeInPath(1); 
  branch firstBranch = constructBranch( modifiedPath.getNthNodeInPath(0), modifiedPath.getNthNodeInPath(2)); 

  modifiedPath.multiplyBranch(traln, rand, firstBranch, multiplier, hastings, prior, brPr); 
  
  /* treat all branches except the first 2 and the last one */
  int ctr = 0; 
  for(nat i = 0; i < modifiedPath.size() ; ++i)
    {
      if(ctr < 2 )
	continue; 

      modifiedPath.multiplyBranch(traln, rand, modifiedPath.at(i), multiplier, hastings, prior, brPr); 
    }

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  modifiedPath.multiplyBranch(traln, rand, constructBranch(sTNode, lastNode), multiplier, hastings, prior, brPr); 
  modifiedPath.multiplyBranch(traln, rand, constructBranch(sTNode, s2LastNode), multiplier, hastings, prior, brPr); 
}



/**
   @brief undoes topological changes, if an eSPR move was done
   according to rPath.

   Only resets the spr move, not any branch length changes due to BL
   multiplying.
 */ 
void SprMove::resetAlongPathForESPR(TreeAln &traln, PriorBelief &prior, const Path& modifiedPath) const
{
  Path anotherPath ; 
  getPathAfterMove(traln, modifiedPath, anotherPath);
  applyPathAsESPR(traln, anotherPath);     
}



/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void SprMove::applyPathAsESPR(TreeAln &traln, const Path &modifiedPath ) const 
{
  tree *tr = traln.getTr();

#ifdef CONTROL_ESPR
double treeLengthBefore = traln.getTreeLengthExpensive(); 
#endif
  
  assert(modifiedPath.size() > 2 ); 

  /* get the subtree ptr */
  nodeptr sTPtr = findNodeFromBranch(tr,getThirdBranch(tr, modifiedPath.at(0), modifiedPath.at(1))); 
  assert(nodeIsInBranch(sTPtr->number, modifiedPath.at(1))); 

  // finds the two nodeptrs adjacent to the subtree  
  nodeptr prOuterPtr = findNodeFromBranch(tr, constructBranch(modifiedPath.getNthNodeInPath(0) ,modifiedPath.getNthNodeInPath(1))),
    prInnerPtr = findNodeFromBranch(tr, constructBranch(modifiedPath.getNthNodeInPath(2), modifiedPath.getNthNodeInPath(1))); 

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  nodeptr s2lPtr = findNodeFromBranch(tr, constructBranch(s2LastNode, lastNode)), // s2l 
    lPtr = findNodeFromBranch(tr, constructBranch(lastNode, s2LastNode)); // l

  nodeptr toBeInserted = prInnerPtr->back;   
  nodeptr toBeInsertedPath = prOuterPtr->back;  
  assert(toBeInserted->number == sTPtr->number && sTPtr->number == toBeInsertedPath->number); 

  /* prune sTPtr */
  traln.clipNode(prInnerPtr, prOuterPtr, prOuterPtr->z[0]); 
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL; 

  traln.clipNode(toBeInsertedPath, lPtr, lPtr->z[0]); 
  traln.clipNode(toBeInserted, s2lPtr, toBeInserted->z[0]); 

  assert(branchExists(tr, constructBranch(sTPtr->number, s2lPtr->number))); 
  assert(branchExists(tr, constructBranch(sTPtr->number, lPtr->number))); 


#ifdef CONTROL_ESPR
  double treeLengthAfter =  traln.getTreeLengthExpensive(); 
  if( fabs(treeLengthAfter  -  treeLengthBefore) > 1e-6  )
    {
      cout << setprecision(8)  << "TL before " << branchLengthToReal(traln.getTr(),treeLengthBefore) << "\tafter" <<  branchLengthToReal(traln.getTr(), treeLengthAfter) << endl; 
      assert(treeLengthAfter == treeLengthBefore); 
    }
#endif
}



/** 
    @brief Gets the description path after the move has been executed. 

    Can be used for reversal => orientation is inverted 

    Does not contain branch lengths.   
 */ 
void SprMove::getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath) const 
{  
  int subTreeNode =  modifiedPath.getNthNodeInPath(1);   
  branch lastBranch = modifiedPath.at(modifiedPath.size()-1); 

  resultPath.append(constructBranch( lastBranch.thatNode , subTreeNode ) )  ; 
  resultPath.append(constructBranch( lastBranch.thisNode , subTreeNode ) )  ; 
  
  assert(modifiedPath.size() > 2); 
  
  // insert the non-descriptive branches in reversed order 
  for(int i = int(modifiedPath.size())-2 ; i > 1 ; --i)
    resultPath.append(modifiedPath.at(i)); 

  // fuse the first two branches 
  auto b1 = modifiedPath.at(0),
    b2  = modifiedPath.at(1); 

  set<int> ids; 
  ids.insert(b1.thisNode); 
  ids.insert(b2.thisNode); 
  ids.insert(b1.thatNode); 
  ids.insert(b2.thatNode); 
  
  assert(ids.find(subTreeNode) != ids.end());   
  ids.erase(subTreeNode); 
  
  assert(ids.size() == 2); 

  
  auto it = ids.begin(); 
  int a = *it; 
  ++it; 
  int b = *it;   
  
  // auto newBranch = constructBranch(a,b); 
  // if(branchesAreConnected(newBranch,  resultPath.at(resultPath.size()-1)

  resultPath.append(constructBranch(a,b)); 

  assert(modifiedPath.size() == resultPath.size()); 
}


ostream& operator<<(ostream &out, const SprMove& rhs)
{
  return out << rhs.path;   
} 
