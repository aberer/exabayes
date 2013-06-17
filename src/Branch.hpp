#ifndef _BRANCH_NEW_H
#define _BRANCH_NEW_H

#include <cassert>
#include "axml.h"
#include "branch.h"

class Branch
{
public: 
  Branch(branch b); 
  Branch(nodeptr p); 
  Branch(nat a = 0, nat b = 0, double length = 0.0); 
  void initFromLegacy(branch b) ; 
  double getInterpretedLength(const TreeAln &traln) const; 
  branch toLegacyBranch() const ; 
  void invert() { swap(thisNode, thatNode) ; }
  Branch getInverted() const { return Branch(thatNode, thisNode, length); }
  bool equalsUndirected(const Branch &rhs) const ;   
  nat getPrimNode() const {return thisNode; } 
  nat getSecNode() const {return thatNode; }
  void setLength(double intLength){length = intLength; }
  double getLength () const {return length; }

  void applyToTree( TreeAln &traln) const ; 

  bool exists(TreeAln &traln) const ;  

  bool nodeIsInBranch(nat node) const {return (thisNode == node) || (thatNode == node) ;  }

  nodeptr findNodeFromBranch(const TreeAln &traln) const; 

  friend ostream& operator<<(ostream &out, const Branch& br); 

  nat getCommonNode(const Branch &rhs ) const; 
  nat getOtherNode(nat node) const {assert(node == thisNode || node ==  thatNode) ; return thisNode == node ? thatNode : thisNode; }

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
    return hash<nat>()(b.getPrimNode()) ^ hash<nat>()( b.getSecNode()) ; 
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
