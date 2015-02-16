#include "Topology.hpp"
#include "common.h"

#include <cassert>
#include <algorithm>

using std::stack; 
using iterator = BareTopology::iterator;

int iterator::ctr = 0; 


iterator::iterator(BareTopology const * ref, Link l) 
  : _ref{ref}
  , _curLink{l}
  , _descent{ }
  , _first{l}
  , _haveBothSides{false}
  , _numTraversed{1} 
{
  
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


void iterator::reset()
{
  auto tmp = stack<Link>{};
  std::swap(_descent, tmp);
}


iterator& iterator::next()
{
  if( _numTraversed == (_ref->size()) )
    *this = _ref->end();
  else
    {
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
                  haveNewOne = true; // actually not, but we terminate here 
                }
              else
                {
                  _curLink = _first.invert();
                  _haveBothSides = true; 
                  haveNewOne = false; 
                }
            }
        }
  
      ++_numTraversed;
    }

  return *this; 
}


iterator& iterator::advance(difference_type n)
{
  if( size_t(_numTraversed) > _ref->size() )
    {
      *this = _ref->end();  
    }
  else if( n > 0 )
    {
      auto t = dynamic_cast< Topology const* >(_ref); 
      if(  t != nullptr )
        {
          auto toJump = size_t(0u);
          // auto &map = t->_bvs;

          if( _curLink.isOuterBranch() )
            {
              toJump = 1u; 
              ++*this;
            }
          else
            {
              toJump = getNumLinksInSubtree();

              if( toJump <= n)
                {
                  if( not _descent.empty() && _curLink == _descent.top())
                    _descent.pop();
                  
                  _curLink = *neighbor();
                  _numTraversed += toJump ;
                  
                  if(  ( _curLink == _first || _curLink == _first.invert() ) 
                       || ( not _descent.empty() &&
                            ( _curLink == _descent.top() || _curLink.invert() == _descent.top() ) ) )
                    {
                      --_numTraversed;
                      ++*this; 
                    }
                }
              else
                {
                  toJump = 1;
                  ++*this;
                }
            }

          *this = advance(n-toJump);
        }
      else
        {
          for(auto i = 0u; i < n; ++i)
            ++*this;
        }
    }

  return *this; 
}


iterator::difference_type distance(iterator first, iterator last)
{
  auto result = iterator::difference_type(0);

  auto t = dynamic_cast<Topology const*>(first._ref);
  if(t != nullptr)
    {
      auto target = t->getBipOrDummy(last);
      auto current = t->getBipOrDummy(first);

      assert(first != t->end() );

      if(  target == current || current == ~target )
        {         
          // => go to return
        }
      else if( first->isOuterBranch() && first->primary() == first->getTaxonNode() )
        {
          // this case avoids skipping the entire tree 
          ++ first;
          result += 1 + distance(first,last); 
        }
      else if( ( (target & current).count() == 0 || ( current <  target)  ) )
        {
          auto n = first.getNumLinksInSubtree();
          first.advance(n);
          result += n + distance(first,last);

        }
      else
        {
          ++first; 
          result += 1 + distance(first, last); 
        }
    }
  else
    {
      while( first != last
             && first != last.opposite() // this is technically incorrect ... 
             ) 
        {
          ++result;
          ++first;
        }
    }

  return result; 
}




size_t iterator::getNumOuterNodesInSubtree() const
{
  auto result = size_t(0); 
  auto t = dynamic_cast<Topology const*>(_ref);
  if( t != nullptr)
    {
      if( _curLink.isOuterBranch() )
        {
          result = _curLink.primary() == _curLink.getTaxonNode()
            ? t->outerSize() - 1
            : 1 ;
        }
      else
        {
          auto found = t->_bvs.find( **this);
          auto foundAlt = t->_bvs.find( (*this)->invert() ) ; 
          if(found  != t->_bvs.end())
            {
              result =  found->second.count();
            }
          else if(  foundAlt != t->_bvs.end() )
            {
              result = t->outerSize() - foundAlt->second.count();
            }
          else
            {
              std::cout << "found neither " << **this << " nor " << (*this)->invert() << " in\n" << t << std::endl;
              assert(0); 
            }
        }
    }
  else
    {
      assert(0);       
    }
  
  return result; 
}



size_t iterator::getNumLinksInSubtree() const
{
  auto result = size_t(0);
  auto t = dynamic_cast<Topology const*>(_ref);
  if( t != nullptr)
    {
      result = getNumOuterNodesInSubtree() * 2 - 1 ; 
    }
  else
    {
      assert(0); 
    }
  
  return result; 
}


bool iterator::isValid() const
{
  auto &map =  _ref->_connections.at(_curLink.primary() ); 
  return map.find(_curLink.secondary()) != map.end() ; 
}


vector<node_id> iterator::findPathTo( iterator end) const
{
  auto result = vector<node_id>{};

  auto casted = dynamic_cast<Topology const*>(_ref); 
  if( casted != nullptr)
    {
      auto endBip = end.getBipOrDummy();
      result = findPathTo_helperWithBipartitions( endBip ) ;
      std::reverse( result.begin(), result.end() );
    }
  else
    {
      result = findPathTo_helper( end );

      if( result.empty() )
        {
          auto tmp = opposite();
          result = tmp.findPathTo_helper(end); 
        }

      assert(not result.empty()); 

      std::reverse(result.begin(), result.end());
  
      result.erase(result.begin() + result.size() - 1  ); 
    }
  
  return result; 
}



vector<node_id> iterator::findPathTo_helperWithBipartitions( bitvector const& endBip) const
{
  auto result = vector<node_id>{};

  auto t = dynamic_cast<Topology const*>(_ref); 
  
  auto curBip = t->getBipOrDummy(*this); 

  auto scoreBip = []( bitvector const& start, bitvector const &end )
    {
      if( start == end)
        return size_t(0) ; 
      else if( start < end)
        return start.symmetricDifference(end).count();
      else if( start < ~end)
        return start.symmetricDifference(~end).count();
      else
        return start.size(); 
    }; 

  auto curScore = scoreBip(curBip, endBip);

  if(curScore == 0 )
    {
      // we are done 
    }
  else
    {
      auto oppo = opposite(); 
      auto oppoScore = scoreBip( ~curBip, endBip);

      auto n1Outer = ( _curLink.isOuterBranch() && _curLink.primary() == _curLink.getTaxonNode() );
      auto n2Outer = ( _curLink.isOuterBranch() && _curLink.primary() == _curLink.getTaxonNode() ); 
      
      auto n1 =  n1Outer ?  opposite().neighbor().opposite()  : neighbor().opposite() ;
      auto n1Bip = t->getBipOrDummy(n1);
      auto n1Score = scoreBip( n1Bip, endBip);

      auto n2 =  n2Outer ? opposite().neighbor().neighbor().opposite() : neighbor().neighbor().opposite();
      auto n2Bip = t->getBipOrDummy(n2);
      auto n2Score = scoreBip(n2Bip,endBip); 

      auto scores = vector<long unsigned>{ oppoScore, n1Score, n2Score };
      auto minScore = std::min_element( scores.begin(), scores.end());

      assert( *minScore < curScore); 
  
      if(  *minScore == oppoScore )
        {
          result = oppo.findPathTo_helperWithBipartitions(endBip);
        }
      else if( *minScore == n1Score )
        {
          result = n1.findPathTo_helperWithBipartitions(endBip);
          result.push_back(n1->secondary()); 
        }
      else if( *minScore == n2Score)
        {
          result = n2.findPathTo_helperWithBipartitions(endBip);
          result.push_back(n2->secondary()); 
        }
      else
        assert(0); 
    }

  return result; 
}




// I do not see a strong reason, why this should be a member function ... 

vector<node_id> iterator::findPathTo_helper( iterator end) const
{
  auto result = vector<node_id>{}; 
  
  if( _curLink == end._curLink || _curLink == end._curLink.invert() )
    {
      result.push_back( _curLink.primary()) ; 
    }
  else if( _curLink.isOuterBranch()  && _curLink.secondary() == _curLink.getTaxonNode())
    {
      // dont add 
    }
  else
    {
      auto n1 = opposite().neighbor();
      auto n2 = n1.neighbor(); 
      
      result = n1.findPathTo_helper(end);

      if( not result.empty() )
        {
          result.push_back(n1->primary());
        }
      else 
        {
          result = n2.findPathTo_helper(end);
          if(not result.empty() )
            result.push_back(n2->primary()); 
        }
    }
  
  return result; 
}


bitvector iterator::getBipOrDummy() const
{
  auto casted = dynamic_cast<Topology const*>(_ref);
  assert(casted != nullptr);
  return casted->getBipOrDummy( *this ); 
}
