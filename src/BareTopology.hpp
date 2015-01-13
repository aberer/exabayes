#ifndef BARE_TOPOLOGY_H
#define BARE_TOPOLOGY_H

#include <map>
#include <stack>
#include <string>
#include <set>
#include <iostream>
#include <vector>
#include <exception>
#include <memory>

#include "bitvector.hpp"

#include "Link.hpp"

class Move; 

// iteratior included below!

class IllegalTreeOperation : public std::exception 
{
public:
  IllegalTreeOperation( std::string str)
    : _str{str}
  {}

  virtual const char* what() const noexcept
  {
    return _str.c_str();
  }
  
  virtual ~IllegalTreeOperation(){}
private:
  std::string _str; 
};


class BareTopology
{
public:
  using node_id = int; 
  
public: 
  class iterator;
  friend class iterator; 
  
  friend std::ostream& operator<<(std::ostream&, const BareTopology&); 
  
public:
  BareTopology();
  /** 
      @brief initializes the tree from a sequence of numbers 
      
      Each number indicates the link, where the next taxon shall be
      inserted. These positions are relative to an iterator that
      starts at a reference taxon (1). The reference taxon allows the
      Topology specialization to rapidly determine these numbers.
  */ 
  BareTopology(  std::vector<size_t> const& numbers);
  BareTopology( std::vector<bitvector> const& bipartitions ); 
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
  iterator end() const;               
  iterator begin( node_id taxon) const;
  /** 
      @brief decomposes the tree and yields a minimal description that
      allows to reconstruct the topology 
      @result a sequence of insertion locations
   */ 
  std::vector<size_t> decomposeToInsertionNumbers(); 
  /**
     @brief inserts the next outer node at a given branch 
   */ 
  virtual iterator insert(iterator it, node_id id = 0 );
  /** 
      @brief removes a subtree represented by an iterator from the
      tree. This subtree must not hold more than one taxon.
      
      @return returns an iterator to the branch existing after removal
   */ 
  virtual iterator erase(iterator it);

  /** 
      @brief applies the move to the tree 
      @return the inverse move 
   */
  virtual std::unique_ptr<Move> move( Move& theMove); 
  
  bool print(std::ostream &s, iterator c, bool needTrifurcation ) const; 
  
  bool operator==(const BareTopology &other) const;

  // for DEBUG 
  void dumpConnections() const ;


  
protected:
  static bool isInnerNode(node_id id )  { return id < 0;  }
  static bool isOuterNode(node_id id) {return id > 0 ; }
  void eraseNode( node_id id); 

  node_id findHighestNode() const  ;
  void initializeWithBipartitions( vector<bitvector> bipartitions );  
  std::tuple<vector<Link>,vector<Link>> getEmergeVansihLinks( iterator movedSubtree, iterator regraftLocation) const;
  
private:
  iterator addBipartitionAt(bitvector curBip, vector<bitvector>  &bipartitions ,iterator branch);
  node_id createInnerNode() { return --_numInnerNodes;   }
  node_id createOuterNode(node_id id = 0 ); 
  
  void hook(Link l);
  void hook( node_id i, node_id o); 
  void unhook(Link l);

private:
  std::map<node_id, std::set<node_id>> _connections;
  int _numInnerNodes;
  int _numOuterNodes;
};



#include "BareTopologyIterator.hpp"


#endif /* TOPOLOGY_H */


