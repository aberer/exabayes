#include "Topology.hpp"

#include "common.h"

#include <algorithm>
#include <cassert>


using iterator = BareTopology::iterator;
using std::unique_ptr;
using std::make_pair;

#include "Move.hpp"


Topology::Topology()
  : BareTopology()
  , _bvs{}
{  
}


std::ostream& operator<<(std::ostream& s, const Topology& c)
{
  auto& t =  static_cast<BareTopology const&>(c); 
  s << t << std::endl;
  for(auto &v : c._bvs)
    s << v.first << "\t" << v.second << std::endl ;
  return s;
}


void Topology::checkedInsert( Link link, bitvector bip)
{
  assert(BareTopology::isOuterNode(link.primary()) || BareTopology::isOuterNode(link.secondary()));
  _bvs.insert(make_pair(link,bip));
}


bitvector Topology::getBipOrDummy(iterator it) const 
{
  auto result = bitvector{};
  
  auto foundA = _bvs.find(*it);
  auto foundB = _bvs.find(it->invert()); 
  
  if( it->isOuterBranch() )
    result.set(it->getTaxonNode() -1 ); 
  else if(  foundA != _bvs.end()  )
    {
      result = foundA->second; 
    }
  else if(  foundB  != _bvs.end() )
    {
      result = ~ foundB->second;
    }
  else
    {
      std::cout << "could not find bipartition for " << *it << std::endl;
      assert(0);
    }
  
  return result; 
}


void Topology::reorient_helper(iterator it , bool forceExistence)
{
  if( it->isOuterBranch()       // is trivial 
      || _bvs.find(*it) != _bvs.end() ) // already there 
    return;
  else
    {
      auto sib =  it.opposite().neighbor();
      reorient_helper( sib, forceExistence );
      reorient_helper( sib.neighbor() , forceExistence);

      if(forceExistence)
        assert( _bvs.find( *(it.opposite()) ) != _bvs.end()  );

      _bvs.erase( *(it.opposite()) );
      
      // correctly insert it 
      auto bip = getBipOrDummy(sib) | getBipOrDummy(sib.neighbor ());
      _bvs.insert(   make_pair(*it , bip) );
    }
}


void Topology::reorient(iterator it, bool forceExistence )
{
  if(outerSize() > 3 )
    {
      if( isOuterNode(it->primary()) )
        it = it.opposite(); 
      
      reorient_helper(it, forceExistence );
      reorient_helper(it.neighbor(), forceExistence);
      reorient_helper(it.neighbor().neighbor(), forceExistence);
    }
}


iterator Topology::erase(iterator it)
{
  auto result = it;

  reorient(it, true);

  auto neighA = it.neighbor();
  auto neighB = neighA.neighbor();

  auto isOuterA = neighA->isOuterBranch();
  auto isOuterB = neighB->isOuterBranch();

  auto bip = bitvector(); 
  
  if(not isOuterA)
    {

      auto foundA = _bvs.find(* neighA);
      bip = foundA->second; 
      assert(foundA != _bvs.end());
      // std::cout << "erasing " << foundA->first << "\t" << foundA->second << std::endl;
      _bvs.erase( foundA);
    }
  
  if(not isOuterB)
    {
      auto foundB = _bvs.find(* neighB);
      bip = foundB->second; 
      assert(foundB != _bvs.end());
      // std::cout << "erasing " << foundB->first << "\t" << foundB->second << std::endl;
      _bvs.erase(foundB);
    }
  
  result = BareTopology::erase(it);

  if( not (isOuterB || isOuterA))
    {
      bip.unset(it->getTaxonNode());
      _bvs.insert( make_pair(*result, bip));
    }
  
  // std::cout << "inserting " << *result << "\t" << bip << std::endl;

  return result;
}



iterator Topology::insert(iterator pos, node_id givenId)
{
  reorient(pos, true);

  auto newPos  = BareTopology::insert(pos, givenId); 

  // reinsert prior link
  if(  outerSize() > 4 && not pos->isOuterBranch() )
  {
    assert(_bvs.find(*pos) != _bvs.end()); 
    auto bip = _bvs.at(*pos);
    
    _bvs.erase(*pos); 
    auto newLink  = Link( newPos->secondary(), pos->secondary() );
    // std::cout << "REINSERT bip in " << newLink << std::endl;
    _bvs[newLink] = bip;
  }

  // insert current link
  if( outerSize() > 3 )
    {
      auto found = newPos.opposite().neighbor();
      auto newOne = found.neighbor();

      if( _bvs.find(*found) == _bvs.end()
          && not found->isOuterBranch() // tagged on; better location possible; 
          )
        std::swap(found, newOne);

      auto bip = getBipOrDummy(found);
      bip.set(newPos->getTaxonNode() -1 );
      _bvs[ *(newOne.opposite())] = bip;
    }
  
  auto highest = findHighestNode() ;
  for(auto &v : _bvs)
    v.second.resize( highest );

  return newPos;
}  


bool Topology::verifyBipartitions() const 
{
  bool okay = true; 
  
  for(auto &v : _bvs)
    {
      auto &link = v.first; 
      auto &bip = v.second;

      auto ctr = bip.count();
      auto iter = iterator(this, link);
      auto controlBv =  bitvector(); 
      while(ctr > 0 )
        {
          if( iter->isOuterBranch() )
            {
              --ctr;
              controlBv.set(iter->getTaxonNode() - 1 ); 
            }
          ++iter; 
        }

      if( controlBv != bip)
        {
          // std::cout << "error: \tcontrol:\t" << controlBv << "\tinthere:\t" << bip  << std::endl;
          okay = false; 
        }
      else
        {
          // std::cout << "GREAT: \tcontrol:\t" << controlBv << "\tinthere:\t" << bip  << std::endl;
        }
    }
  return okay; 
}





bool Topology::isEquivalent( Topology const& rhs) const
{
  // this requires some computational effort, but given the data
  // structure this is equivalence test that is easiest to implement.
  
  auto cpyA = *this;
  auto cpyB = rhs;

  cpyA.reorient( cpyA.begin(1), true);
  cpyB.reorient(cpyB.begin(1), true ); 
  
  auto setUnion = vector<bitvector>();
  auto setInter = vector<bitvector>();
  
  auto setA = std::set<bitvector>();
  auto setB = std::set<bitvector>();

  auto lam = []( std::pair<Link,bitvector> const& elem ){ return elem.second; }; 
  std::transform(cpyA._bvs.begin(), cpyA._bvs.end(), std::inserter(setA, setA.begin()),lam);
  std::transform(cpyB._bvs.begin(), cpyB._bvs.end(), std::inserter(setB, setB.begin()),lam);

  std::set_union( setA.begin(), setA.end(), setB.begin(), setB.end(), std::back_inserter(setUnion) );
  std::set_intersection( setA.begin(), setA.end(), setB.begin(), setB.end(), std::back_inserter(setInter));

  return setInter.size() == setUnion.size() ; 
}


unique_ptr<Move> Topology::move(Move &theMove )
{
  auto van = theMove.getVanishingLinks();
  
  auto movedSubtree = theMove.getReferenceBranch(); 
  
  reorient(movedSubtree, true );
  
  BareTopology::move(theMove);
  
  for(auto v : van )
    {
      auto iter = _bvs.find( v);
      auto iter2 = _bvs.find( v.invert()); 
      if(  iter  != _bvs.end() )
        _bvs.erase( iter);
      else if( iter2 != _bvs.end() )
        _bvs.erase(iter2); 
    }

  reorient(movedSubtree, false);

  return theMove.getInverse(); 
}
