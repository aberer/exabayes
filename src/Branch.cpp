#include "Branch.hpp"
#include "TreeAln.hpp"


Branch::Branch(nat a , nat b, double length) 
  : thisNode(a), thatNode(b), length(length)
{
}


bool Branch::isTipBranch(const TreeAln &traln) const
{ 
  nat  num = traln.getTr()->mxtips ; 
  return (thisNode <= num) || (thatNode <= num);
}


bool Branch::exists(const TreeAln &traln) const
{
  nodeptr p = traln.getTr()->nodep[thisNode]; 
  return
    (nat)p->back->number == thatNode
    || (nat)p->next->back->number == thatNode
    || (nat)p->next->next->back->number == thatNode; 
}



nodeptr Branch::findNodePtr(const TreeAln &traln) const
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



std::ostream& operator<<(std::ostream &out, const Branch& br)
{ 
  return out << "(" << br.thisNode << "/" << br.thatNode << "):" <<  std::setprecision(std::numeric_limits<double>::digits10 ) << br.length; 
}



bool Branch::equalsUndirected(const Branch &rhs) const 
{ 
  return ( rhs.getPrimNode() == thisNode && rhs.getSecNode() == thatNode  ) 
    ||  ( rhs.getPrimNode( )== thatNode && rhs.getSecNode() == thisNode ) ; 
}



// TODO replace as well 
double Branch::getInterpretedLength(const TreeAln &traln) const
{ 
  assert(traln.getNumBranches() == 1 ); 
  return -log(length) * traln.getTr()->fracchange;  
} 


double Branch::getInternalLength(const TreeAln &traln, double length) const
{
  assert(traln.getNumBranches() == 1 ); 
  return exp( - length / traln.getTr()->fracchange) ; 
}


// nat Branch::getCommonNode(const Branch &rhs ) const
// {
//   if(thisNode == rhs.thisNode || thatNode == rhs.thisNode)
//     return rhs.thisNode;
//   else if(thatNode == rhs.thatNode || thisNode == rhs.thatNode)
//     return rhs.thatNode; 
//   else 
//     return 0;   
// }



Branch Branch::getThirdBranch(const TreeAln &traln, const Branch& rhs ) const
{
  // TODO efficiency 
  int node = getIntersectingNode(rhs); 
  assert(not traln.isTipNode(traln.getTr()->nodep[node])); 

  nodeptr p = traln.getNode(node); 
  nodeptr q = p->back,
    q1 = p->next->back,
    q2 = p->next->next->back; 
  
  int altA = q->number,
    altB = q1->number,
    altC = q2->number; 
  
  int pNumber = p->number; 
  
  if( not rhs.nodeIsInBranch(altA) &&   not nodeIsInBranch(altA ) )    
    return Branch(pNumber, altA) ; 
  else if( not   rhs.nodeIsInBranch(altB) &&   not nodeIsInBranch(altB) )
    return Branch(pNumber,altB) ;
  else if( not   rhs.nodeIsInBranch(altC) &&   not nodeIsInBranch(altC) )
    return Branch(pNumber, altC) ;
  else 
    {
      assert(0); 
      return Branch(0,0); 
    }
} 


nat Branch::getIntersectingNode(const Branch  &rhs) const 
{
  if(rhs.nodeIsInBranch( thisNode ) )
    return thisNode; 
  else if(rhs.nodeIsInBranch(thatNode))
    return thatNode; 
  else 
    {
      assert(0); 
      return 0; 
    }
} 


void Branch::readFromCheckpoint( std::istream &in )
{
  thisNode = cRead<nat>(in); 
  thatNode = cRead<nat>(in);   
  length = cRead<double>(in); 
} 

void Branch::writeToCheckpoint( std::ostream &out)  const
{
  cWrite(out, thisNode); 
  cWrite(out, thatNode); 
  cWrite(out, length); 
}  

