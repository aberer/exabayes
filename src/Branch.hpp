/** 
    @file Branch.hpp 

    @brief represents a branch 

    Furthermore, Branch may also represent (depending on context) a
    sub-tree or be equivalent to a nodeptr. In these cases, a
    secondary node defines the orientation of the branch. 
 */ 

#ifndef _BRANCH_NEW_H
#define _BRANCH_NEW_H

#include <iostream>
#include <cassert>
#include <vector>

#include "axml.h"
#include "Checkpointable.hpp"

class AbstractParameter; 
class TreeAln; 


enum class BranchEqualFlag : int 
{
  WITHOUT_ANYTHING = 0,
  WITH_LENGTH = 1, 
    WITH_DIRECTION = 2
}; 


BranchEqualFlag operator|( BranchEqualFlag a, BranchEqualFlag b); 
BranchEqualFlag operator&( BranchEqualFlag a, BranchEqualFlag b); 


class Branch : public Checkpointable
{
public: 
  /////////////////
  // life cycle  //
  /////////////////
  Branch(nat a = 0, nat b = 0, std::vector<double> lengths =std::vector<double>{}); 

  ///////////////
  // observers //
  ///////////////
  /** 
      @brief gets the absolute (true) length of the branch
   */ 
  double getInterpretedLength(const TreeAln &traln, const AbstractParameter* param) const; 

  // double getInternalLength(const TreeAln &traln, double length) const; 
  void setConvertedInternalLength(const TreeAln& traln, const AbstractParameter* param, double length) ; 
  /**
     @brief switches the orientation of the branch (i.e., the reference node )
   */ 
  Branch getInverted() const { return Branch(thatNode, thisNode, lengths); }
  /**
     @brief indicates whether two branches are equal disregarding branch lengths and orientation     
   */ 
  bool equalsUndirected(const Branch &rhs) const ;   
  /** 
      @brief gets the equality of this branch and rhs in a generic manner. 
   */ 
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
      @brief gets all lengths stored in the branch
   */ 
  const std::vector<double>& getAllLengths() const {return lengths; }
  /**
     @brief sets the primary node
   */ 
  void setPrimNode(nat node)  {thisNode = node; }
  /** 
      @brief sets the secondary (reference) node
   */ 
  void setSecNode(nat node) {thatNode = node; }  
  /**
     @brief sets the branch length (internal representation)
   */ 
  void setLength(double intLength, const AbstractParameter* param);
  /**
     @brief gets the (internal) branch length
   */ 
  double getLength (const AbstractParameter* param) const;
  /**
     @brief indicates whether the branch exists in tree traln
   */ 
  bool exists(const TreeAln &traln) const; 
  /**
     @brief gets the branch adjacent to the branch "this" and "rhs"
   */ 
  Branch getThirdBranch(const TreeAln &traln, const Branch& rhs ) const; 
  /**
     @brief gets the intersecting node of two branches
     @notice triggers assertion, if there is no intersection 
  */ 
  nat getIntersectingNode(const Branch  &rhs) const ; 
  /** 
      @brief indicates whether a node is part of this branch 
   */ 
  bool nodeIsInBranch(nat node) const {return (thisNode == node) || (thatNode == node) ;  }
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

  virtual void readFromCheckpoint( std::istream &in ) ; 
  virtual void writeToCheckpoint( std::ostream &out) const  ;   

  // friends 
  friend std::ostream& operator<<(std::ostream &out, const Branch& br); 


  static bool printLength; 

private: 
  nat thisNode; 
  nat thatNode;   
  std::vector<double> lengths; 
}; 




class BranchHashNoLength
{ 
public: 
  size_t operator()(const Branch &b ) const 
  {
    return std::hash<nat>()(b.getPrimNode()) ^ std::hash<nat>()( b.getSecNode()) ; 
  }
}; 

class BranchEqualNoLength
{
public:
  bool operator() (const Branch &b, const Branch &a) const 
  {
    return a.equalsUndirected(b); 
  }
};

#endif
