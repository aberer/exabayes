#include "Branch.hpp"
#include "TreeAln.hpp"
#include "parameters/AbstractParameter.hpp"


Branch::Branch(nat a , nat b, std::vector<double> _lengths) 
  : thisNode(a), thatNode(b), lengths(_lengths)
{
}


BranchEqualFlag operator|( BranchEqualFlag a, BranchEqualFlag b) 
{
  return static_cast<BranchEqualFlag>(static_cast<int>(a) | static_cast<int>(b)); 
}


BranchEqualFlag operator&( BranchEqualFlag a, BranchEqualFlag b)
{
  return static_cast<BranchEqualFlag>(static_cast<int>(a) & static_cast<int>(b)); 
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



std::ostream& operator<<(std::ostream &out, const Branch& rhs)
{ 
  out << MAX_SCI_PRECISION; 
  out << "(" << rhs.thisNode << "/" << rhs.thatNode
	     << "):[" ; 
  for(auto &v : rhs.lengths)
    out << v << ","; 
  out << "]";
  return out; 
}



bool Branch::equals(const Branch &rhs, BranchEqualFlag flags) const 
{
  bool result = true; 

  // dir a 
  result &= (rhs.thisNode == thisNode)  && (rhs.thatNode == thatNode); 
  
  if((flags & BranchEqualFlag::WITH_DIRECTION) == BranchEqualFlag::WITHOUT_ANYTHING )
    result |= (rhs.thatNode == thisNode) && (rhs.thisNode == thatNode); 

  if( (flags & BranchEqualFlag::WITH_LENGTH) != BranchEqualFlag::WITHOUT_ANYTHING)
    {
      for(nat i = 0; i < rhs.lengths.size(); ++i)
	result &= (rhs.lengths[i] == lengths[i]) ;  
    }

  return result; 
}



bool Branch::equalsUndirected(const Branch &rhs) const 
{ 
  return equals(rhs, BranchEqualFlag::WITHOUT_ANYTHING); 
}


double Branch::getInterpretedLength(const TreeAln &traln, const AbstractParameter* param) const
{ 
  double fracC = traln.getMeanSubstitutionRate(param->getPartitions()); 

  // TODO that's potentially dangerous 
  double length = 0; 
  if(lengths.size() == 1 )
    length = lengths.at(0);
  else 
    length = lengths.at(param->getIdOfMyKind()); 
  
  return -log(length)* fracC; 
} 


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
  
  nat num = cRead<nat>(in);
  for(nat i = 0; i < num ; ++i)
    lengths.push_back(cRead<double>(in));
} 

void Branch::writeToCheckpoint( std::ostream &out)  const
{
  cWrite(out, thisNode); 
  cWrite(out, thatNode); 
  cWrite(out, lengths.size());
  for(auto &v : lengths)
    cWrite(out, v); 
}  


void Branch::setConvertedInternalLength(const TreeAln& traln, 
					const AbstractParameter* param, double length) 
{
  double fracC = traln.getMeanSubstitutionRate(param->getPartitions());  
  double internalLength = exp(- length / fracC); 
  setLength(internalLength, param);
} 

void Branch::setLength(double intLength, const AbstractParameter* param)
{
  lengths.resize(MAX(lengths.size(), param->getIdOfMyKind() + 1 ), TreeAln::initBL); 
  lengths.at(param->getIdOfMyKind()) = intLength; 
}


double Branch::getLength (const AbstractParameter* param) const
{
  return lengths.at(param->getIdOfMyKind());
}
