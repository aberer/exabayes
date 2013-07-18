#ifndef _BRANCH_NEW_H
#define _BRANCH_NEW_H

#include <iostream>
#include <cassert>

#include "axml.h"
#include "Checkpointable.hpp"

class TreeAln; 

class Branch : public Checkpointable
{
public: 
  /////////////////
  // life cycle  //
  /////////////////
  Branch(nat a = 0, nat b = 0, double length = 0.0); 

  ///////////////
  // observers //
  ///////////////
  /** 
      @brief gets the absolute (true) length of the branch
   */ 
  double getInterpretedLength(const TreeAln &traln) const; 
  /**
     @brief converts the absolute length "length" into the internal
     respresentation and sets it
   */
  double getInternalLength(const TreeAln &traln, double length) const; 
  /**
     @brief switches the orientation of the branch (i.e., the reference node )
   */ 
  Branch getInverted() const { return Branch(thatNode, thisNode, length); }
  /**
     @brief indicates whether two branches are equal disregarding branch lengths and orientation     
   */ 
  bool equalsUndirected(const Branch &rhs) const ;   
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
     @brief sets the branch length (internal representation)
   */ 
  void setLength(double intLength){length = intLength; }
  /**
     @brief gets the (internal) branch length
   */ 
  double getLength () const {return length; }
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
  
  virtual void readFromCheckpoint( std::ifstream &in ) ; 
  virtual void writeToCheckpoint( std::ofstream &out)  ;   

  // friends 
  friend std::ostream& operator<<(std::ostream &out, const Branch& br); 

private: 
  nat thisNode; 
  nat thatNode;   
  double length; 
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
