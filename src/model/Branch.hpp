/** 
    @file Branch.hpp 

    @brief represents a branch 

    Branch comes in three variants: without branch length
    (BranchPlain), with a single branch length (BranchLength, usefull
    when e.g. operating on a single partition) and all branch lengths
    possible for this branch (BranchLengths, usefull when backing up
    data)

    The current implementation relies on some hard-coding for the
    three instances -- this could be improved by using more templates.

 */ 

#ifndef _BRANCH_NEW_NEW_H
#define _BRANCH_NEW_NEW_H

#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <unordered_set>

#include <type_traits>

#include "pll.h"

#include "system/Serializable.hpp"

class TreeAln; 
#include "system/FlagType.hpp" 

class AbstractParameter; 

enum class BranchEqualFlag : int 
{
  WITHOUT_ANYTHING = 0,
  WITH_LENGTH = 1, 
    WITH_DIRECTION = 2
}; 


template<typename TYPE>
class LengthPart; 

template<typename TYPE = void>
class Branch : public Serializable, public LengthPart<TYPE>
{
public: 
  /////////////////
  // life cycle  //
  /////////////////
  Branch(nat a = 0, nat b = 0)
    : thisNode(a)
    , thatNode(b)
  {
  };

  ///////////////
  // observers //
  ///////////////
  /**
     @brief switches the orientation of the branch (i.e., the reference node )
   */ 
  Branch getInverted() const ; 
  /**
     @brief indicates whether two branches are equal disregarding branch lengths and orientation     
   */ 
  bool equalsUndirected(const Branch &rhs) const ;   
  /** 
      @brief gets the equality of this branch and rhs in a generic manner. 
   */ 
  // TODO replace this with a proper operator== once in a while ... 
  bool equals(const Branch &rhs, BranchEqualFlag flags) const ;  
  /**
     @brief gets the primary node (i.e., reference node) of the branch
   */ 
  nat getPrimNode() const {return thisNode; } 
  /**
     @brief gets the secondary node (that is there for orientation purposes)
   */ 
  nat getSecNode() const {return thatNode; }  
  /**
     @brief sets the primary node
   */ 
  void setPrimNode(nat node)  {thisNode = node; }
  /** 
      @brief sets the secondary (reference) node
   */ 
  void setSecNode(nat node) {thatNode = node; }  
  /**
     @brief indicates whether the branch exists in tree traln
   */ 
  bool exists(const TreeAln &traln) const; 
  /**
     @brief gets the branch adjacent to the branch "this" and "rhs"
   */ 
  Branch<void> getThirdBranch(const TreeAln &traln, const Branch& rhs ) const; 
  /**
     @brief gets the intersecting node of two branches
     @notice triggers assertion, if there is no intersection 
  */ 
  nat getIntersectingNode(const Branch  &rhs) const ; 
  /** 
      @brief indicates whether a node is part of this branch 
  */ 
  bool hasNode(nat node) const {return (thisNode == node) || (thatNode == node) ;  }
  /** 
      @brief finds the nodeptr in the tree that is described by this
      branch
   */ 
  nodeptr findNodePtr(const TreeAln &traln) const; 
  /** 
      @brief gets the other node (given a node)
   */ 
  nat getOtherNode(nat node) const {assert(node == thisNode || node ==  thatNode) ; return thisNode == node ? thatNode : thisNode; }
  /** 
      @brief indicates whether a branch is a tip branch  
   */ 
  bool isTipBranch(const TreeAln &traln) const; 
  /**
     @brief gets the distance to another branch
   */ 
  nat getDistance(const Branch &bOther, const TreeAln &traln) const ; 
  nat getDistanceHelper(const Branch &bOther, const TreeAln &traln, nat distSoFar ) const ; 
  /** 
      @brief gets the taxon ids on one side of a branch (which equals one partition of a bipartition)
   */ 
  std::vector<nat> getBipartition(const TreeAln &traln ) const ; 

  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const  ;   


  Branch<double> toBlDummy() const; 
  Branch<std::vector<double>> toBlsDummy() const ; 
  Branch<void> toPlain() const ; 

  Branch<double> toOneLength(const AbstractParameter* param) const; 

  friend std::ostream& operator<< (std::ostream &out, const Branch<TYPE>& rhs)
  {
    out << "(" << rhs.getPrimNode() << "," << rhs.getSecNode() << ")"; 
    rhs.lengthToString(out); 
    return out; 
  }
 
  bool isAdjacent(const Branch<void> &rhs) const ;  

private: 			// METHODS 
  std::unordered_set<nat> getBipartitionHelper(const TreeAln &traln  ) const ; 


private: 			// ATTRIBUTES 
  nat thisNode; 
  nat thatNode;   
  std::vector<double> lengths; 
};

#include "LengthPart.hpp"


typedef Branch<void> BranchPlain; 
typedef Branch<double> BranchLength; 
typedef Branch<std::vector<double>> BranchLengths; 


namespace std 
{
  template<> struct hash<Branch<void>>
  {
    size_t operator()(const Branch<void> & x) const
    {
      return std::hash<size_t>()(x.getSecNode()) ^ std::hash<size_t>()(x.getPrimNode()); 
    }
  };

  template<> struct hash<Branch<double>>
  {
    size_t operator()(const Branch<double> & x) const
    {
      return hash<Branch<void>>()(x.toPlain());
    }
  }; 

  template<> struct hash<Branch<std::vector<double>>>
  {
    size_t operator()(const Branch<std::vector<double>> & x) const
    {
      return hash<Branch<void>>()(x.toPlain());
    }
  }; 

  template<> struct equal_to<Branch<void>>
  {
    bool operator()(const Branch<void> &a, const Branch<void> &b)  const
    {
      return a.equalsUndirected(b); 
    }
  }; 

  template<>
  struct equal_to<Branch<double>>
  {
    bool operator()(const Branch<double> &a, const Branch<double> &b)const
    {
      return a.equalsUndirected(b); 
    }
  }; 

  template<>
  struct equal_to<Branch<std::vector<double>>>
  {
    bool operator()(const Branch<std::vector<double>> &a, const Branch<std::vector<double>> &b)const
    {
      return a.equalsUndirected(b); 
    }
  }; 
}

std::ostream& operator<< (std::ostream &out, const Branch<void>& rhs); 
std::ostream& operator<< (std::ostream &out, const Branch<double>& rhs); 
std::ostream& operator<< (std::ostream &out, const Branch<std::vector<double>>& rhs); 

#endif
