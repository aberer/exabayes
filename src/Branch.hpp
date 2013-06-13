#ifndef _BRANCH_NEW_H
#define _BRANCH_NEW_H

#include "axml.h"
#include "branch.h"


class Branch
{
public: 
  Branch(nat a = 0, nat b = 0, double length = 0.0); 
  void initFromLegacy(branch b) ; 
  double getInterpretedLength(const TreeAln &traln) const; 
  branch toLegacyBranch() const ; 
  void invert() { swap(thisNode, thatNode) ; }
  bool equalsUndirected(const Branch &rhs) const ; 
  nat getPrimNode() const {return thisNode; } 
  nat getSecNode() const {return thatNode; }

  friend ostream& operator<<(ostream &out, const Branch& br); 

private: 
  nat thisNode; 
  nat thatNode; 
  
  double length; 

}; 




#endif
