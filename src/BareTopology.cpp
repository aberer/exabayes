#include "BareTopology.hpp"

#include <iostream>
#include <cassert>
#include "NotImplementedException.hpp"

using std::make_pair;

using iterator = BareTopology::iterator;

BareTopology::BareTopology()
  : _connections{}
  , _numInnerNodes{0}
  , _numOuterNodes{0}
{
}


void BareTopology::unhook(Link l)
{
  // std::cout <<  "UNHOOK " << l  << std::endl;
  
  auto p = l.primary();
  auto s = l.secondary(); 
  
  auto &mapP =  _connections[p];
  auto &mapS = _connections[s]; 
  
  auto iter = mapP.find(s);
  assert(iter != mapP.end());
  mapP.erase(iter);
  
  iter = mapS.find(p);
  assert(iter != mapS.end());
  mapS.erase(iter);
}


void BareTopology::hook(node_id i, node_id o)
{
  assert( outerSize() == 2 || not (isOuterNode(i) && isOuterNode(o)));
  _connections[i].insert(o);
  _connections[o].insert(i);
}


iterator BareTopology::insert(iterator it)
{
  auto result = end();
  if( outerSize() == 0 )
    {
      auto id = createOuterNode();
      _connections[id] = set<node_id>{}; 
      result = iterator(this, Link(id, 0)); 
    }
  else if ( outerSize() == 1 )
    {
      auto oId = createOuterNode();
      auto o = it._curLink._priNode; 
      hook(oId, o);
      result = iterator(this, Link(oId, o)); 
    }
  else
    {
      auto iId = createInnerNode();
      auto oId = createOuterNode();
      hook(oId, iId);
      hook( iId, it._curLink._priNode );
      hook( iId, it._curLink._secNode );
      result = iterator(this, Link(oId, iId)); 
    }

  // if *it is fully initialized
  if(it->primary() != 0 && it->secondary() != 0 )
    unhook(*it);

  return result;
}


iterator BareTopology::begin( node_id taxon ) const
{
  if( outerSize() == 1  )
    return iterator(this, Link(taxon,0)); 
  else
    {
      auto other = _connections.at(taxon).begin();
      return iterator(this , Link(taxon, *other));
    }
}


iterator BareTopology::begin() const 
{
  auto result = end(); 

  if(outerSize() == 0 )
    result = end(); 
  else if(outerSize() == 1 )
    {
      result = iterator(this, Link(1,0)); 
    }
  else 
    {
      auto elem = _connections.begin();
      auto beg = elem->first ; 
      auto other =
        // elem->second.begin() == elem->second.end()
        // ? 0 :
      *(elem->second.begin()); 
      result =  iterator(this, Link(beg,other)); 
    }
  
  return result;  
}


iterator BareTopology::end() const
{
  return iterator(this, Link{}); 
}


iterator iterator::neighbor() const 
{
  if( _curLink.primary() > 0 )
    return *this;

  auto result = *this; 
  auto& map = result._ref->_connections.at(result._curLink._priNode); 
  auto found = map.find(result._curLink._secNode);
  assert(found != map.end());
  ++found;

  // cycle around 
  if(found == map.end())
    found = map.begin();
  
  result._curLink._secNode = *found;
  result.reset();
  return result; 
}


iterator iterator::opposite( ) const 
{
  auto result = *this;  
  result._curLink = result._curLink.invert();
  result.reset();
  return result; 
}


iterator::iterator(BareTopology const * ref, Link l) 
  : _ref{ref}
  , _curLink{l}
  , _descent{ }
  , _first{l}
  , _haveBothSides{false}
{
  
}



iterator& iterator::next()
{
  // std::cout << "NEXT of " << **this << std::endl;
  
  auto haveNewOne = false;
  while(not haveNewOne)
    {
      _curLink = isOuterNode( _curLink.secondary())
        ? *neighbor()
        : *opposite().neighbor(); 
        
      haveNewOne = _descent.empty() || ( not (_curLink == _descent.top().invert() || _curLink == _descent.top())) ;
      if(not haveNewOne)
        _descent.pop();
      else if(not _curLink.isOuterBranch())
        _descent.push(_curLink);

      if( _curLink == _first.invert() || _curLink == _first )
        {
          if( _haveBothSides  || _first.isOuterBranch() )
            {
              *this = _ref->end();
              break;
            }
          else
            {
              _curLink = _first.invert();
              _haveBothSides = true; 
              haveNewOne = false; 
            }
        }
    }
  return *this; 
}


iterator& iterator::advance(size_t n)
{
  // TODO short-cut? 
  
  for(auto i = 0u; i < n; ++i)
    ++*this;
  return *this; 
}



std::ostream& operator<<(std::ostream& s, const Link& c)
{
  s << c._priNode << "," << c._secNode; 
  return s;
}



bool BareTopology::print(std::ostream &s, iterator c, bool needTrifurcation ) const
{
  if( c->primary() > 0 )
    {
      s << c->getTaxonNode();
      return true; 
    }
  else
    {
      auto sib1 = c.neighbor().opposite();
      auto sib2  = c.neighbor().neighbor().opposite();

      if(not needTrifurcation)
        s << "("  ;
      print(s,sib1, false); 
      s<<   "," ;
      print(s, sib2, false) ;
      if(not needTrifurcation)
        s << ")";
      return false; 
    }
}


std::ostream& operator<<(std::ostream& s, BareTopology const& c)
{
  if(c.outerSize() == 0)
    {
      s << "();"; 
    }
  else
    {
      auto start = c.begin();
      auto needTrifurcation = true; 
      s << "(" ;
      needTrifurcation = c.print(s,start, needTrifurcation) ;
      s << "," ;
      c.print(s,start.opposite(), needTrifurcation ) ;
      s << ");";
    }

  return s;
}



void iterator::reset()
{
  auto tmp = stack<Link>{};
  std::swap(_descent, tmp);
}
  
