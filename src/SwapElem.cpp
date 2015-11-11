#include "SwapElem.hpp"

#include "comm/ParallelSetup.hpp"

nat SwapElem::getMyId(ParallelSetup& pl, nat runid) const 
{
  assert(not pl.isMyChain(runid, _idA) != not pl.isMyChain(runid, _idB)); 
  if(pl.isMyChain(runid, _idA))
    return _idA; 
  else 
    return _idB; 
}

nat SwapElem::getRemoteId(ParallelSetup& pl, nat runid) const 
{
  nat mine = getMyId(pl,runid); 
  if(mine == _idA)
    return _idB; 
  else 
    return _idA; 
}

bool operator==(const SwapElem &elemA, const SwapElem &elemB ) 
{
  bool isEqual = true; 
  isEqual &=    (elemA._idA ==  elemB._idA && elemA._idB ==  elemB._idB)
    || (elemA._idA ==  elemB._idB && elemA._idA ==  elemB._idB); 
  isEqual &= elemA._gen == elemB._gen; 
  isEqual &= elemA._r == elemB._r; 
  return isEqual; 
}
  
std::ostream& operator<<(std::ostream& out, const SwapElem &elem)
{
  out << "gen=" << elem._gen << " " <<  elem._idA << "<->"  <<  elem._idB  <<  " " << elem._r; 
  return out; 
}



bool SwapElem::onlyOneIsMine(ParallelSetup& pl, nat runid) const 
{
  auto oneIsMine = ( not pl.isMyChain(runid, _idA) ) != ( not  pl.isMyChain(runid,_idB)) ; 
  return oneIsMine; 
}
 


bool SwapElem::bothAreMine(ParallelSetup& pl, nat runid) const 
{
  auto bothAreMine = pl.isMyChain(runid, _idA) && pl.isMyChain(runid,_idB); 
  return bothAreMine; 
}
