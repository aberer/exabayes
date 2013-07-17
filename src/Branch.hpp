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
  double getInterpretedLength(const TreeAln &traln) const; 
  double getInternalLength(const TreeAln &traln, double length) const; 
  Branch getInverted() const { return Branch(thatNode, thisNode, length); }
  bool equalsUndirected(const Branch &rhs) const ;   
  nat getPrimNode() const {return thisNode; } 
  nat getSecNode() const {return thatNode; }  
  void setPrimNode(nat node)  {thisNode = node; }
  void setSecNode(nat node) {thatNode = node; }  
  void setLength(double intLength){length = intLength; }
  double getLength () const {return length; }
  bool exists(const TreeAln &traln) const; 
  Branch getThirdBranch(const TreeAln &traln, const Branch& rhs ) const; 
  nat getIntersectingNode(const Branch  &rhs) const ; 
  bool nodeIsInBranch(nat node) const {return (thisNode == node) || (thatNode == node) ;  }
  nodeptr findNodePtr(const TreeAln &traln) const; 
  nat getCommonNode(const Branch &rhs ) const; 
  nat getOtherNode(nat node) const {assert(node == thisNode || node ==  thatNode) ; return thisNode == node ? thatNode : thisNode; }
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
