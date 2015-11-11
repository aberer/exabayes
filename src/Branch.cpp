#include "Branch.hpp"
#include "TreeAln.hpp"
#include "parameters/AbstractParameter.hpp"
#include <limits>


template<>
void Branch<std::vector<double> >::deserialize( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   

  assert(0); 
  // TODO lengths 
}

template<>
void Branch<double>::deserialize( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   
  length =  cRead<decltype(length)>(in);
} 


template<>
void Branch<void>::deserialize( std::istream &in )
{
  thisNode = cRead<decltype(thisNode)>(in); 
  thatNode = cRead<decltype(thatNode)>(in);   
}

template<>
void Branch<double>::serialize( std::ostream &out)  const
{
  cWrite<decltype(thisNode)>(out, thisNode); 
  cWrite<decltype(thatNode)>(out, thatNode); 
  cWrite<decltype(length)>(out, length); 
}

template<>
void Branch<std::vector<double>>::serialize( std::ostream &out)  const
{
  assert(0); 
}


template<>
void Branch<void>::serialize( std::ostream &out)  const
{
  cWrite<decltype(thisNode)>(out, thisNode); 
  cWrite<decltype(thatNode)>(out, thatNode); 
}  



template<>
Branch<void> Branch<void>::getInverted() const
{
  auto result = Branch<void>(thatNode,thisNode); 
  return result; 
}

template<>
Branch<double> Branch<double>::getInverted() const
{
  auto result = Branch<double>(thatNode,thisNode); 
  result.length = this->length; 
  
  // if(std::is_same<BT, std::vector<double>>::value)
  //   result.lengths = this->lengths; 

  return result; 
} 

template<>
Branch<std::vector<double>> Branch<std::vector<double>>::getInverted() const
{
  auto result = Branch<std::vector<double>>(thatNode,thisNode); 
  result.lengths = this->lengths; 
  return result; 
}



template<>
Branch<double> Branch<std::vector<double>>::toOneLength(const AbstractParameter* param) const
{
  auto result = Branch<double>(thisNode,thatNode); 
  result.setLength(getLength(param)); 
  return result; 
} 





template<typename TYPE>
bool Branch<TYPE>::isTipBranch(const TreeAln &traln) const
{ 
  nat  num = traln.getTrHandle().mxtips ; 
  return (thisNode <= num) || (thatNode <= num);
}


template<typename TYPE>
bool Branch<TYPE>::exists(const TreeAln &traln) const
{
  nodeptr p = traln.getTrHandle().nodep[thisNode]; 
  return
    (nat)p->back->number == thatNode
    || (nat)p->next->back->number == thatNode
    || (nat)p->next->next->back->number == thatNode; 
}

template<typename TYPE>
nodeptr Branch<TYPE>::findNodePtr(const TreeAln &traln) const
{
  auto& tr = traln.getTrHandle(); 
  
  nodeptr p = tr.nodep[thisNode]; 
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


template<typename TYPE>
bool Branch<TYPE>::equals(const Branch &rhs, BranchEqualFlag flags) const 
{
  bool result = true; 

  // dir a 
  result &= (rhs.thisNode == thisNode)  && (rhs.thatNode == thatNode); 
  
  if((flags & BranchEqualFlag::WITH_DIRECTION) == BranchEqualFlag::WITHOUT_ANYTHING )
    result |= (rhs.thatNode == thisNode) && (rhs.thisNode == thatNode); 

  
  if(std::is_same<TYPE, std::vector<double>>::value || std::is_same<TYPE,double>::value)
    {
      if( (flags & BranchEqualFlag::WITH_LENGTH) != BranchEqualFlag::WITHOUT_ANYTHING)
	{
	  for(nat i = 0; i < rhs.lengths.size(); ++i)
	    result &= (rhs.lengths[i] == lengths[i]) ;  
	}
    }

  return result; 
}


template<typename TYPE>
bool Branch<TYPE>::equalsUndirected(const Branch &rhs) const 
{ 
  return equals(rhs, BranchEqualFlag::WITHOUT_ANYTHING); 
}


template<typename TYPE>
BranchPlain Branch<TYPE>::getThirdBranch(const TreeAln &traln, const Branch& rhs ) const
{
  int node = getIntersectingNode(rhs); 
  assert(not traln.isTipNode(traln.getTrHandle().nodep[node])); 

  nodeptr p = traln.getNode(node); 
  nodeptr q = p->back,
    q1 = p->next->back,
    q2 = p->next->next->back; 
  
  int altA = q->number,
    altB = q1->number,
    altC = q2->number; 
  
  int pNumber = p->number; 
  
  if( not rhs.hasNode(altA) &&   not hasNode(altA ) )    
    return BranchPlain(pNumber, altA) ; 
  else if( not   rhs.hasNode(altB) &&   not hasNode(altB) )
    return BranchPlain(pNumber,altB) ;
  else if( not   rhs.hasNode(altC) &&   not hasNode(altC) )
    return BranchPlain(pNumber, altC) ;
  else 
    {
      assert(0); 
      return BranchPlain(0,0); 
    }
} 


template<typename TYPE>
nat Branch<TYPE>::getIntersectingNode(const Branch  &rhs) const 
{
  if(rhs.hasNode( thisNode ) )
    return thisNode; 
  else if(rhs.hasNode(thatNode))
    return thatNode; 
  else 
    {
      assert(0); 
      return 0; 
    }
} 



template<typename TYPE>
auto Branch<TYPE>::getBipartition(const TreeAln &traln ) const -> std::vector<nat>
{
  auto bip = this->getBipartitionHelper(traln);
  auto complement = this->getInverted().getBipartitionHelper(traln); 

  auto result = std::vector<nat>(); 


  if(isTipBranch(traln))
    {
      result.insert(result.end(), bip.begin(), bip.end()); 

    }
  else 
    { 
      assert(bip.size() + complement.size() == traln.getNumberOfTaxa()); 

      if(bip.size( )< complement.size() )
	result.insert(result.end(), bip.begin(), bip.end()); 
      else 
	result.insert(result.end(), complement.begin(), complement.end()); 
    }
  
  std::sort(result.begin(), result.end(), std::less<nat>()); 

  return result; 
}



template<typename TYPE>
auto Branch<TYPE>::getBipartitionHelper(const TreeAln &traln  ) const -> std::unordered_set<nat> 
{
  auto result = std::unordered_set<nat>();

  if(  isTipBranch(traln) )
    {
      if(traln.isTipNode(findNodePtr(traln)))
	result.insert(getPrimNode()); 
      else 
	result.insert(getSecNode()); 
    }
  else 
    {

#if 0 
      // just commented out for compilation speed 
      auto desc = traln.getDescendents(*this); 
      auto bipA = desc.first.getInverted().getBipartitionHelper(traln);
      auto bipB = desc.second.getInverted().getBipartitionHelper(traln);
      result.insert( bipA.begin(), bipA.end()); 
      result.insert( bipB.begin(), bipB.end()); 

#else
      assert(0); 
#endif

    }

  return result; 
}


template<typename TYPE>
nat Branch<TYPE>::getDistanceHelper(const Branch &bOther, const TreeAln &traln, nat distSoFar ) const 
{
  // tout << "checking " << toString() << " against " << bOther.toString() << std::endl; 

  if(equalsUndirected(bOther)) 
    {
      // tout << "found and dist is " << distSoFar << std::endl; 
      return distSoFar; 
    }
  else if ( isTipBranch(traln) && traln.isTipNode(findNodePtr(traln)))
    {
      // tout << toString() << " is tip branch, did notfind " << std::endl; 
      return std::numeric_limits<nat>::max();
    }
  else 
    {
      auto desc = traln.getDescendents(this->toPlain()); 
      auto distA = desc.first.getInverted().getDistanceHelper(bOther.toPlain(), traln, distSoFar+1); 
      auto distB = desc.second.getInverted().getDistanceHelper(bOther.toPlain(), traln, distSoFar + 1 ); 
      return std::min(distA, distB); 
    }
}

template<typename TYPE>
nat Branch<TYPE>::getDistance(const Branch &bOther, const TreeAln &traln ) const 
{
  auto distA = getDistanceHelper(bOther, traln,0 ); 
  auto distB = getInverted().getDistanceHelper(bOther, traln, 0); 
  auto trueDist = std::min(distA,distB); 
  assert(trueDist < std::numeric_limits<nat>::max()); 
  assert(std::max(distA,distB) == std::numeric_limits<nat>::max() 
	 || ( distA == 0 && distB == 0)  );
  return trueDist; 
}


template<typename T>
BranchPlain Branch<T>::toPlain() const  
{
  return BranchPlain(this->thisNode, this->thatNode); 
}

template<typename T>
Branch<double> Branch<T>::toBlDummy() const
{
  return Branch<double>{thisNode, thatNode}; 
}



template<typename T>
Branch<std::vector<double>> Branch<T>::toBlsDummy() const 
{
  return Branch<std::vector<double>> {thisNode, thatNode}; 
}

  

template<typename T>
bool Branch<T>::isAdjacent(const BranchPlain &rhs) const 
{
  return rhs.hasNode(thisNode) || rhs.hasNode(thatNode); 
}


template class Branch<void>;
template class Branch<double>;
template class Branch<std::vector<double>>;


  
