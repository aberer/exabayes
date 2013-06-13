#include "Branch.hpp"
#include "TreeAln.hpp"



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
