#include "BareTopology.hpp"
#include "Topology.hpp"

#include "common.h"

#include <sstream>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "all_exceptions.hpp"

#include "Move.hpp"


using std::tuple; 
using std::string; 
using std::vector; 
using std::stack;
using std::map; 
using std::set; 
using std::make_pair;

using iterator = BareTopology::iterator;

#include "SPRMove.hpp"

BareTopology::BareTopology()
  : _connections{}
  , _numInnerNodes{0}
  , _numOuterNodes{0}
{
}


BareTopology::BareTopology(  std::vector<size_t> const& numbers)
  : BareTopology()
{
  // default initialization 
  insert(begin());
  insert(begin());
  insert(begin());

  for(auto &v : numbers)
    {
      assert(v <= size()); 
      insert(begin(1) + v ); 
    }
}



BareTopology::BareTopology( vector<bitvector> const& bipartitions )
  : BareTopology()
{
  initializeWithBipartitions(bipartitions);
}



// TODO no copy 
void BareTopology::initializeWithBipartitions( vector<bitvector>  bipartitions )
{
  // TODO use a reverse iterator instead 

  std::reverse(bipartitions.begin(), bipartitions.end()); 
  // we must assume that bipartitions have been created via a
  // post-order travarsal. The last bipartition (now the first) *must*
  // contain all taxa.

  assert(bipartitions.front().count() == bipartitions.front().size() ); 

  // auto innerNode = begin(1)->secondary(); 

  // addFirstBipartitionAt( bipartitions, innerNode);
  auto topBip = bipartitions.front();
  bipartitions.erase(bipartitions.begin()); 
  addBipartitionAt(topBip, bipartitions, begin()); 
}




iterator BareTopology::addBipartitionAt(bitvector curBip, vector<bitvector>&  bipartitions, iterator branch)
{
  // std::cout << "ADD "<< curBip << " to "  << *branch  << std::endl;
  auto alreadyAdded = bitvector();
  alreadyAdded.resize(curBip.size()); // NECESSARY??? 
  
  // treat sub-bipartitions
  auto branchForNext =  branch; 
  while( not bipartitions.empty() && bipartitions.front() < curBip )
    {
      auto nextBip = bipartitions.front();
      alreadyAdded |= nextBip;
      bipartitions.erase(bipartitions.begin());
      branchForNext = addBipartitionAt(nextBip, bipartitions, branchForNext); 
    }

  
  if( alreadyAdded < curBip)
    {
      auto toInsert = vector<size_t> (curBip - alreadyAdded);
      assert(toInsert.size() == 1 );

      // this branch must exist
      // std::cout << "INSERT\t" << toInsert[0] << "\t" << *branch  << std::endl;
      insert(branch, static_cast<node_id>(toInsert[0] + 1 ));
    }

  
  // return the correct position of the bipartition in the tree so far
  auto p  = branch->primary();
  auto s = branch->secondary(); 

  if(  p ==  0 || s == 0 || size() <  2 )
    {
      branch = begin();
    }
  else if(  _connections.at(p).find(s) == _connections.at(p).end() ) 
    {
      assert( _connections.find(p) != _connections.end()); 
      auto &map = _connections.at(p);
      assert( _connections.find(s) != _connections.end()); 
      auto &oppoMap = _connections.at(s);
      auto res = vector<node_id>();

      std::set_intersection( map.begin(), map.end(),  oppoMap.begin(), oppoMap.end(), std::back_inserter(res));
      assert(res.size() == 1 );

      auto &conOfRes = _connections.at(res[0]);
      auto other = vector<node_id>();
      auto known = std::set<node_id>{p, s }; 
      std::set_difference(conOfRes.begin(), conOfRes.end(),   known.begin(), known.end(), std::back_inserter(other) ) ; 
      assert(other.size() == 1 );

      branch._curLink = Link(res[0] , other[0]);
    }
  else
    {
      std::cout << "no modification necessary" << std::endl; 
    }

  return branch; 
}


void BareTopology::unhook(Link l)
{
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



void BareTopology::hook( node_id i, node_id o)
{
  hook(Link(i,o));
}



void BareTopology::hook(Link l)
{
  node_id i = l.primary(); 
  node_id o = l.secondary(); 

  assert( outerSize() == 2 || not (isOuterNode(i) && isOuterNode(o)));
  _connections[i].insert(o);
  _connections[o].insert(i);
}





iterator BareTopology::insert(iterator it, node_id givenId  )
{
  auto result = end();
  if( outerSize() == 0 )
    {
      auto id = createOuterNode( givenId );
      _connections[id] = set<node_id>{}; 
      result = iterator(this, Link(id, 0)); 
    }
  else if ( outerSize() == 1 )
    {
      auto oId = createOuterNode( givenId);
      auto o = it._curLink._priNode; 
      hook( Link(oId, o) );
      result = iterator(this, Link(oId, o)); 
    }
  else
    {
      auto iId = createInnerNode( );
      auto oId = createOuterNode(givenId);
      hook(Link(oId, iId));
      hook( Link(iId, it._curLink._priNode) );
      hook( Link(iId, it._curLink._secNode )  );
      result = iterator(this, Link(oId, iId)); 
    }

  // if *it is fully initialized
  if(it->primary() != 0 && it->secondary() != 0 )
    unhook(*it);

  return result;
}


void BareTopology::eraseNode( node_id id)
{
  if( id <  0)
    ++_numInnerNodes;
  else
    --_numOuterNodes;

  assert( _connections[id].empty() );
  _connections.erase( _connections.find( id )); 
}


iterator BareTopology::erase(iterator it)
{
  auto result = it;
  if( outerSize() == 1 )
    {
      eraseNode(1);
    }
  else if( size() == 1 )
    {
      assert( Link(1,2) == *it); 
      unhook( Link(1,2));
      eraseNode(2);
      result = iterator(this, Link(1,0)); 
    }
  else if( not (
                it->isOuterBranch()
                &&  it->secondary() == it->getTaxonNode()
                && it->secondary() == static_cast<int>(outerSize())
                )
           )
    {
      throw UnsupportedTreeAction();
    }
  else
    {
      auto itA = it.neighbor();
      auto itB = itA.neighbor();
  
      unhook(*it);
      unhook(*itA);
      unhook(*itB);
      
      auto l = Link(itA->secondary(), itB->secondary()); 
  
      eraseNode(it->primary());
      eraseNode(it->secondary());

      hook(   l.primary() , l.secondary()   );
  
      result.reset();
      result._curLink = l ; 
    }

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
      result = iterator(this, Link(findHighestNode(),0)); 
    }
  else 
    {
      auto elem = _connections.begin();
      auto beg = elem->first ; 
      auto other = *(elem->second.begin()); 
      result =  iterator(this, Link(beg,other)); 
    }
  
  return result;  
}


iterator BareTopology::end() const
{
  return iterator(this, Link{}); 
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
  else if( c.outerSize() == 1 )
    {
      s << "1;"  ; 
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


bool BareTopology::operator==(const BareTopology &other) const
{
  for(auto &val : _connections)
    {
      auto &f = val.first;
      auto &s = val.second;

      if(  other._connections.find(f) == other._connections.end() )
        return false; 

      auto &s2  = other._connections.at(f);

      if( s.size() != s2.size())
        return false;

      for(auto &v : s)
        {
          if( s2.find(v) == s2.end())
            return false; 
        }
    }
  
  return true; 
}



void BareTopology::dumpConnections() const  
{
  for(auto iter = _connections.begin(); iter != _connections.end(); ++iter)
    {
      auto &from = iter->first ;
      std::cout << from << "\t->\t";
      for(auto i = iter->second.begin(); i != iter->second.end(); ++i)
        {
          std::cout << *i <<",";
        }
      std::cout << "\n";
    }
}



vector<size_t> BareTopology::decomposeToInsertionNumbers()
{
  auto result = vector<size_t>();

  while( size() > 3 )
    {
      auto highest = findHighestNode();
      
      auto last = begin( highest  ).opposite();
      auto newIt = erase(last);
      auto dist = distance(begin(1), newIt);
      result.push_back( dist );
    }

  std::reverse(result.begin(), result.end()); 
  
  return result; 
}


node_id BareTopology::createOuterNode(node_id id )
{
  auto result = id;       
  if(  id == 0 )
    {
      // TODO  INEFFICIENT 
      
      // determine the id yourself

      bool hasBroken = false; 
      // naive expensive operation 
      for(auto i = node_id(1); i < _numOuterNodes + 1 ; ++i)
        {
          if( _connections.find(i) == _connections.end() )
            {
              hasBroken = true; 
              result = i;
              break; 
            }
        }

      if(not hasBroken)
        result = _numOuterNodes + 1 ;

      assert( _connections.find(result)== _connections.end() && result != 0); 
    }
  else if( _connections.find(id) !=   _connections.end() )
    {
      auto &&ss = std::stringstream{};
      ss << id; 
      throw IllegalTreeOperation("trying to insert outer node " + ss.str()  + " that is already there." ) ; 
    }

  ++_numOuterNodes;
  // std::cout << "created " << result << std::endl;
  return result; 
}


node_id BareTopology::findHighestNode() const
{
  auto const& elem = std::max_element( _connections.begin() , _connections.end() ,
                                       []( std::pair<node_id, std::set<node_id>>  const &elemA ,std::pair<node_id, std::set<node_id>>  const &elemB ) { return elemA.first < elemB.first; });
  return elem->first; 
}


std::unique_ptr<Move> BareTopology::move(Move &theMove  )
{
  for(auto &v : theMove.getVanishingLinks())
    unhook( v );

  for(auto &v : theMove.getEmergingLinks())
    hook(v);

  theMove.apply();

  return theMove.getInverse();
}



