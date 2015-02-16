#ifndef LINK_H
#define LINK_H


#include <iosfwd>

#include <map>

using node_id = int;


class Link
{
  friend class BareTopology;
  friend class Topology; 

public:
  Link( node_id a = 0, node_id b = 0 )
    : _priNode{a}
    , _secNode{b}
  {}

  bool operator==(  Link const& rhs) const { return _priNode == rhs._priNode && _secNode == rhs._secNode; }
  bool operator!=(const Link &other) const { return !(*this == other); }
  
  friend std::ostream& operator<<(std::ostream& s, const Link& c); 

  Link invert() const { return Link(_secNode,_priNode); }

  bool isOuterBranch() const { return _priNode > 0 || _secNode >  0;  }

  node_id getTaxonNode() const { return _priNode < 0 ? _secNode : _priNode;       } 
  
  node_id primary() const {  return _priNode; }
  node_id secondary() const  { return _secNode ; }
  
private: 
  node_id _priNode;             // primary node represented by this link 
  node_id _secNode;             // secondary node represented by this link
};

namespace std
{
  template<>
  struct hash<Link>
  {
    size_t operator()(const Link & x) const
    {
      return std::hash<size_t>()( x.primary()  ) ^ std::hash<size_t>()( x.secondary() ); 
    }
  }; 

  template<> 
  struct equal_to<Link>
  {
    bool operator()(const Link &a, const Link &b)  const
    {
      return a == b ; 
    }
  }; 
}



#endif /* LINK_H */
