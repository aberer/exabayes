#ifndef BARE_TOPOLOGY_H
#define BARE_TOPOLOGY_H

#include <map>
#include <stack>
#include <set>
#include <iostream>

using std::stack;
using std::map; 
using std::set; 

using node_id = int;


class BareTopology;

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



class BareTopology
{
public: 
  class iterator;
  friend class iterator; 
  
  friend std::ostream& operator<<(std::ostream&, const BareTopology&); 
  
public:
  BareTopology();
  virtual ~BareTopology(){}
  /** 
      @brief gets the number of outer nodes 
   */ 
  size_t outerSize() const { return _numOuterNodes; }
  /**
     @brief gets the number of inner nodes 
   */ 
  size_t innerSize() const { return - _numInnerNodes; }
  /** 
      @brief gets the number of branches in the tree 
   */ 
  size_t size() const { return outerSize() + innerSize() -1 ; }
  
  iterator begin() const; 
  iterator end() const;               // grml also const

  iterator begin( node_id taxon) const;

  /**
     @brief inserts the next outer node at a given branch 
   */ 
  virtual iterator insert(iterator it);

  bool print(std::ostream &s, iterator c, bool needTrifurcation ) const; 
  
  void dumpConnections() const 
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
  
protected:
  static bool isInnerNode(node_id id )  { return id < 0;  }
  static bool isOuterNode(node_id id) {return id > 0 ; }

private: 
  node_id createInnerNode() { return --_numInnerNodes;   }
  node_id createOuterNode() { return ++_numOuterNodes; }
  
  void hook(node_id i, node_id o);
  void unhook(Link l);

private:
  map<node_id, set<node_id>> _connections;
  int _numInnerNodes;
  int _numOuterNodes;
};


class BareTopology::iterator
{
  friend class BareTopology;
  friend class Topology;

  using iterator_category = std::forward_iterator_tag ;
  using value_type = Link;
  using difference_type = size_t;
  using pointer = Link*;
  using reference = Link&; 
  
public:
  iterator(BareTopology const*  ref, value_type l);
  iterator(  const iterator& rhs) = default;
  iterator& operator=(const iterator &rhs) = default; 
  /** 
      @brief equals ->next and looses any traversal memory   
   */ 
  iterator neighbor(  ) const;
  /** 
      @brief equals ->back and looses any traversal memory   
  */ 
  iterator opposite( ) const;
  
  iterator& next();

  bool operator==(const iterator &other) const { return _ref == other._ref && _curLink == other._curLink; }
  bool operator!=(const iterator &other) const { return !(*this == other); }
  
  Link const& operator*() const {  return _curLink; }
  Link const *  operator->() const { return &_curLink ;  }
  
  Link get() const {return _curLink; }

  iterator& operator++() { next(); return *this; }
  iterator operator+( size_t n ) { auto cpy = *this; cpy.advance(n); return cpy; }
  iterator& advance(size_t n); 

  /** 
      @brief forgets the branches it already has traversed 
   */ 
  void reset() ;
private:                        // METHODS 
  /** 
      @brief conducts extra steps necessary, once the first subtree is traversed  
  */
  void handleHalfIsDone(bool &haveNewOne) ; 
private:
  BareTopology const* _ref; 
  value_type _curLink;
  stack<value_type> _descent;
  Link _first;
  bool _haveBothSides;
};


#endif /* TOPOLOGY_H */
