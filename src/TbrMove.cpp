#include "TbrMove.hpp"


void TbrMove::applyToTree(TreeAln &traln) const 
{
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1); 

  applyPath(traln,path1); 
  applyPath(traln,path2); 
} 


void TbrMove::revertTree(TreeAln &traln, PriorBelief &prior) const 
{
  Path path1r ; 
  Path path2r ; 
  getPathAfterMove(traln,path1, path1r);
  getPathAfterMove(traln,path2, path2r);
  
  applyPath(traln,path1r);
  applyPath(traln,path2r);

  path1.restoreBranchLengthsPath(traln, prior); 
  path2.restoreBranchLengthsPath(traln, prior); 
} 


void TbrMove::disorientAtNode(TreeAln &traln, nodeptr p) const 
{
  // limitation: must be the bisecting branch currently 

  int nodeA = path1.getNthNodeInPath(1),
    nodeB = path2.getNthNodeInPath(1); 

  nodeptr  q = nullptr, r = nullptr ; 
  if(nodeA == p->number)
    {
      assert(nodeB == p->back->number); 
      q = p; 
      r = p->back; 
    }
  else if( nodeB  == p->number)
    {
      assert(nodeA == p->back->number); 
      q = p->back; 
      r = p ; 
    }
  else 
    {
      assert(0);
    }

  sprDisorientPath(traln,q, path1);
  sprDisorientPath(traln,r, path2);
} 


void TbrMove::extractMoveInfo(const TreeAln &traln, std::vector<Branch> description) 
{
  sprCreatePath(traln, description.at(0),description.at(1),path1); 
  sprCreatePath(traln, description.at(2),description.at(3),path2); 
} 


void TbrMove::multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<AbstractPrior*> brPr) const 
{
  assert(0);
} 


Branch TbrMove::getEvalBranch(const TreeAln &traln) const 
{   
  return  Branch(path1.getNthNodeInPath(1) , path2.getNthNodeInPath(1)); 
}


std::ostream& operator<<(std::ostream &out, const TbrMove &rhs) 
{
  return out <<  "path1:"  << rhs.path1 << std::endl
	     << "path2:" << rhs.path2 << std::endl; 
}  
