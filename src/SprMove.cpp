#include <set>

#include "SprMove.hpp"
#include "Branch.hpp"

// TODO make branch stuff more optional (performance)

// TODO constructor instead of extract move info 


void SprMove::applyToTree(TreeAln &traln) const
{
  applyPath(traln, path); 
}

void SprMove::revertTree(TreeAln &traln, PriorBelief &prior) const
{
  Path anotherPath ; 
  getPathAfterMove(traln, path, anotherPath);
  applyPath(traln, anotherPath); 
  path.restoreBranchLengthsPath(traln, prior); 
}

 
void SprMove::disorientAtNode(TreeAln &traln, nodeptr p) const
{
  sprDisorientPath(traln,p, path);
}


void SprMove::extractMoveInfo(const TreeAln &traln, std::vector<Branch> description)
{
  sprCreatePath(traln, description.at(0), description.at(1), path);
} 


AbstractMove* SprMove::clone() const
{
  return new SprMove; 
}


Branch SprMove::getEvalBranch(const TreeAln &traln) const
{    
  Branch b = path.at(0).getThirdBranch(traln, path.at(1)); 
  Branch lastBranch = path.at(path.size()-1); 

  Branch bA = Branch(b.getPrimNode(), lastBranch.getPrimNode()), 
    bB = Branch(b.getPrimNode(),lastBranch.getSecNode()); 

  assert(bA.exists(traln) && bB.exists(traln)); 

  Branch futureRoot = bA.getThirdBranch(traln, bB ); 
  
  return futureRoot; 
}



/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void SprMove::multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<AbstractPrior*> brPrs)  const 
{  
  assert(path.size() >= 2); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
  assert(brPrs.size() < 2 ); 
  auto brPr = brPrs[0]; 
 
  int sTNode = path.getNthNodeInPath(1); 
  Branch firstBranch = Branch( path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 

  path.multiplyBranch(traln, rand, firstBranch, multiplier, hastings, prior, brPr); 
  
  /* treat all branches except the first 2 and the last one */
  int ctr = 0; 
  for(nat i = 0; i < path.size() ; ++i)
    {
      if(ctr < 2 )
	continue; 

      path.multiplyBranch(traln, rand, path.at(i), multiplier, hastings, prior, brPr); 
    }

  int lastNode = path.getNthNodeInPath(path.getNumberOfNodes()-1),
    s2LastNode = path.getNthNodeInPath(path.getNumberOfNodes()-2); 
  path.multiplyBranch(traln, rand, Branch(sTNode, lastNode), multiplier, hastings, prior, brPr); 
  path.multiplyBranch(traln, rand, Branch(sTNode, s2LastNode), multiplier, hastings, prior, brPr); 
} 


/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void SprMove::applyPath(TreeAln &traln, const Path &modifiedPath ) const 
{
  assert(modifiedPath.size() > 2 ); 

  /* get the subtree ptr */
  
  Branch third =  modifiedPath.at(0).getThirdBranch(traln, modifiedPath.at(1)); 
  nodeptr sTPtr = third.findNodePtr(traln); 
  assert(modifiedPath.at(1).nodeIsInBranch(sTPtr->number )); 

  

  // finds the two nodeptrs adjacent to the subtree  
  nodeptr prOuterPtr = Branch(modifiedPath.getNthNodeInPath(0) ,modifiedPath.getNthNodeInPath(1)).findNodePtr(traln ),
    prInnerPtr = Branch(modifiedPath.getNthNodeInPath(2), modifiedPath.getNthNodeInPath(1)).findNodePtr(traln ); 

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  nodeptr s2lPtr = Branch(s2LastNode, lastNode).findNodePtr(traln ), // s2l 
    lPtr = Branch(lastNode, s2LastNode).findNodePtr(traln ); // l

  nodeptr toBeInserted = prInnerPtr->back;   
  nodeptr toBeInsertedPath = prOuterPtr->back;  
  assert(toBeInserted->number == sTPtr->number && sTPtr->number == toBeInsertedPath->number); 

  /* prune sTPtr */
  traln.clipNode(prInnerPtr, prOuterPtr, prOuterPtr->z[0]); 
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL; 

  traln.clipNode(toBeInsertedPath, lPtr, lPtr->z[0]); 
  traln.clipNode(toBeInserted, s2lPtr, toBeInserted->z[0]); 

  assert(Branch(sTPtr->number, s2lPtr->number).exists(traln)); 
  assert(Branch(sTPtr->number, lPtr->number).exists(traln )); 
}



/** 
    @brief Gets the description path after the move has been executed. 

    Can be used for reversal => orientation is inverted 

    Does not contain branch lengths.   
 */ 
void SprMove::getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath) const 
{  
  int subTreeNode =  modifiedPath.getNthNodeInPath(1);   
  Branch lastBranch = modifiedPath.at(modifiedPath.size()-1); 

  resultPath.append(Branch( lastBranch.getSecNode() , subTreeNode ) )  ; 
  resultPath.append(Branch( lastBranch.getPrimNode() , subTreeNode ) )  ; 
  
  assert(modifiedPath.size() > 2); 
  
  // insert the non-descriptive branches in reversed order 
  for(int i = int(modifiedPath.size())-2 ; i > 1 ; --i)
    resultPath.append(modifiedPath.at(i)); 

  // fuse the first two branches 
  auto b1 = modifiedPath.at(0),
    b2  = modifiedPath.at(1); 

  std::set<int> ids; 
  ids.insert(b1.getPrimNode()); 
  ids.insert(b2.getPrimNode()); 
  ids.insert(b1.getSecNode()); 
  ids.insert(b2.getSecNode()); 
  
  assert(ids.find(subTreeNode) != ids.end());   
  ids.erase(subTreeNode); 
  
  assert(ids.size() == 2); 

  
  auto it = ids.begin(); 
  int a = *it; 
  ++it; 
  int b = *it;   
  
  // auto newBranch = Branch(a,b); 
  // if(branchesAreConnected(newBranch,  resultPath.at(resultPath.size()-1)

  resultPath.append(Branch(a,b)); 

  assert(modifiedPath.size() == resultPath.size()); 
}


std::ostream& operator<<(std::ostream &out, const SprMove& rhs)
{
  return out << rhs.path;   
} 




void SprMove::sprCreatePath(const TreeAln &traln, Branch mover, Branch movedInto, Path &pathHere ) const
{
#ifdef EFFICIENT
  // looks very clunky 
  assert(0); 
#endif

  Branch chosen = movedInto; // description.at(1); 
  Branch prunedTree = mover; // description.at(0); 

  nodeptr p = prunedTree.findNodePtr(traln);
  
  pathHere.clear();   
  pathHere.findPath(traln, p, chosen.findNodePtr(traln));
  pathHere.reverse();   

  Branch Tmp = pathHere.at(0); 

  nodeptr aNode = Tmp.findNodePtr(traln); 
  if(traln.isTipNode(aNode))
    aNode = Tmp.getInverted().findNodePtr(traln);
  assert(not traln.isTipNode(aNode)); 
  Branch b = pathHere.at(0); 
  while(b.equalsUndirected( pathHere.at(0)) || b.equalsUndirected(Branch(p->number, p->back->number))) 
    {
      aNode = aNode->next; 
      b = Branch(aNode->number, aNode->back->number) ;
    }
  pathHere.reverse();
  pathHere.append(b);  
  pathHere.reverse();

  // TODO inefficient
  pathHere.saveBranchLengthsPath(traln); 
}


void SprMove::sprDisorientPath(TreeAln &traln, nodeptr p, const Path &pathHere) const 
{  
  assert(pathHere.size() > 2) ; 
  
  int first = (int)pathHere.getNthNodeInPath(0),
    last = (int)pathHere.getNthNodeInPath(pathHere.getNumberOfNodes()-1); 

  if(not pathHere.nodeIsOnPath(p->number) || traln.isTipNode(p)
     || p->number == first || p->number == last )
    return; 

  disorientHelper(traln, p);
  
  sprDisorientPath( traln, p->next->back, pathHere); 
  sprDisorientPath( traln, p->next->next->back, pathHere);
}
