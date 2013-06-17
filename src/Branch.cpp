#include "Branch.hpp"
#include "TreeAln.hpp"



Branch::Branch(nodeptr p)
  : Branch(p->number, p->back->number, p->z[0])
{
}


Branch::Branch(branch b)
  : Branch(b.thisNode, b.thatNode)
{
  this->length = b.length[0]; 
}


Branch::Branch(nat a , nat b, double length) 
  : thisNode(a), thatNode(b), length(length)
{
}


void Branch::initFromLegacy(branch b) 
{
  this->thisNode = b.thisNode; 
  this->thatNode = b.thatNode; 
  length = b.length[0]; 
}
  

bool Branch::exists(TreeAln &traln) const
{
  // nodeptr p = findNodeFromBranch(traln); 
  nodeptr p = traln.getTr()->nodep[thisNode]; 
  return
    (nat)p->back->number == thatNode
    || (nat)p->next->back->number == thatNode
    || (nat)p->next->next->back->number == thatNode; 
}



nodeptr Branch::findNodeFromBranch(const TreeAln &traln) const
{
  auto tr = traln.getTr(); 
  
  nodeptr p = tr->nodep[thisNode]; 
  if(p->back->number == (int)thatNode)
    return p ; 
  else if(p->next->back->number == (int)thatNode) 
    return p->next; 
  else 
    {
      assert(p->next->next->back->number == (int)thatNode); 
      return p->next->next; 
    }  
} 



ostream& operator<<(ostream &out, const Branch& br)
{ 
  return out << "(" << br.thisNode << "/" << br.thatNode << "):" << br.length; 
}



bool Branch::equalsUndirected(const Branch &rhs) const 
{ 
  return ( rhs.getPrimNode() == thisNode && rhs.getSecNode() == thatNode  ) 
    ||  ( rhs.getPrimNode( )== thatNode && rhs.getSecNode() == thisNode ) ; 
}


double Branch::getInterpretedLength(const TreeAln &traln) const
{ 
  return -log(length) * traln.getTr()->fracchange;  
} 



branch Branch::toLegacyBranch() const
{
  auto b = constructBranch(thisNode,thatNode) ; 
  b.length[0] = length; 
  return b;
}


nat Branch::getCommonNode(const Branch &rhs ) const
{
  if(thisNode == rhs.thisNode || thatNode == rhs.thisNode)
    return rhs.thisNode;
  else if(thatNode == rhs.thatNode || thisNode == rhs.thatNode)
    return rhs.thatNode; 
  else 
    return 0;   
} 



void Branch::applyToTree( TreeAln &traln) const
{
  nodeptr p = findNodeFromBranch(traln); 
  double tmp = length; 
  // =/ 
  traln.clipNode(p, p->back, tmp); 
}
