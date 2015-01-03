#include "Topology.hpp"

#include <cassert>

using std::make_pair; 


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
  if( it->isOuterBranch() )
    result.set(it->getTaxonNode() -1 ); 
  else
    {
      auto found = _bvs.find(*it); 
      result = found->second; 
    }
  
  return result; 
}


void Topology::reorient_helper(iterator it )
{
  if( it->isOuterBranch()       // is trivial 
      || _bvs.find(*it) != _bvs.end() ) // already there 
    return;
  else
    {
      auto sib =  it.opposite().neighbor();
      reorient_helper( sib );
      reorient_helper( sib.neighbor() );

      // delete incorrectly oriented bip 
      assert( _bvs.find( *(it.opposite()) ) != _bvs.end()  );
      // std::cout << "ERASE " << *(it.opposite()) << "\t"  << _bvs[*(it.opposite())] << std::endl;
      _bvs.erase( *(it.opposite()) );

      // correctly insert it 
      auto bip = getBipOrDummy(sib) | getBipOrDummy(sib.neighbor ());
      _bvs.insert(   make_pair(*it , bip) );
      // std::cout << "RE-INSERT " << *it << "\t" << bip << std::endl;
      
    }
}


void Topology::reorient(iterator it)
{
  if(outerSize() > 3 )
    {
      if( isOuterNode(it->primary()) )
        it = it.opposite(); 
      
      reorient_helper(it );
      reorient_helper(it.neighbor());
      reorient_helper(it.neighbor().neighbor());
    }
}


iterator Topology::insert(iterator pos)
{
  reorient(pos);

  auto newPos  = BareTopology::insert(pos); 

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

      if( _bvs.find(*found) == _bvs.end() )
        std::swap(found, newOne);

      auto bip = getBipOrDummy(found);
      bip.set(newPos->getTaxonNode() -1 );
      _bvs[ *(newOne.opposite())] = bip;
      // std::cout << "INSERT new bip in " << *newOne << std::endl;
    }

  // TODO efficiency problem? 
  for(auto &v : _bvs)
    v.second.resize(outerSize() + 1 );

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
