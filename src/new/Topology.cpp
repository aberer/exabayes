#include "Topology.hpp"

#include "LikeArrayManager.hpp"
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
  , _observingLnlArray{}
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
  insertBip(link, bip);
}


bitvector Topology::getBipOrDummy(iterator it) const 
{
  auto result = bitvector{};
  
  auto foundA = _bvs.find(*it);
  auto foundB = _bvs.find(it->invert()); 
  
  if( it->isOuterBranch() )
    {
      result.resize( outerSize() ); 
      result.set(it->getTaxonNode() -1 );
    }
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

      eraseBip(*(it.opposite()) ) ; 

      // correctly insert it 
      auto bip = getBipOrDummy(sib) | getBipOrDummy(sib.neighbor ());

      insertBip( *it, bip);
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
      
      if(_observingLnlArray != nullptr)
        _observingLnlArray->tellRoot(it); 
    }
}


iterator Topology::erase(iterator it)
{
  auto result = it;

  reorient(it, true);

  auto neighA = it.neighbor();
  auto neighB = neighA.neighbor();

  auto bip = bitvector(); 
  
  auto removeIfInner = [&](iterator anIt )
    {
      auto res = not anIt->isOuterBranch(); 
      if( res )
        {
          auto found = _bvs.find( *anIt);
          bip = found->second;
          assert(found != _bvs.end());
          eraseBip(found->first); 
        }
      return res; 
    }; 

  auto isOuterA = removeIfInner(neighA);
  auto isOuterB = removeIfInner(neighB); 
  
  result = BareTopology::erase(it);

  if( not (isOuterB || isOuterA))
    {
      bip.unset(it->getTaxonNode());
      insertBip(*result, bip); 
    }

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

    eraseBip(*pos); 

    insertBip(Link( newPos->secondary(), pos->secondary() ),
              bip); 
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

      insertBip(*(newOne.opposite()),bip);  
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
        eraseBip( iter->first); 
      else if( iter2 != _bvs.end() )
        eraseBip( iter2->first); 
    }

  reorient(movedSubtree, false);

  return theMove.getInverse(); 
}


void Topology::eraseBip( Link const& l )
{
  auto found = _bvs.find(l); 

  if(found != _bvs.end())
    {
      auto bip =  found->second; 
      
      _bvs.erase(found);
  
      if(_observingLnlArray != nullptr)
        _observingLnlArray->makeVanish(bip);
    }
}


void Topology::insertBip( Link l, bitvector b)
{
  _bvs.insert( make_pair(l , b) );
  
  if(_observingLnlArray != nullptr)
  _observingLnlArray->makeExist(b); 
}


void Topology::setObservingLnlArrayManager( LikeArrayManagerPtr ptr)
{
  _observingLnlArray = ptr;
  for(auto &v : _bvs)
    _observingLnlArray->makeExist(v.second);
}
