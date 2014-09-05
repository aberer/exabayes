#ifndef NODEAGE_H
#define NODEAGE_H

#include "BranchPlain.hpp"

// for this class, _primNode is the node, we identify the age for 

// _secNode points upwards (with the NOT having any upward node)


class NodeAge : public BranchPlain
{
public:
  NodeAge(const BranchPlain &b = BranchPlain(), double h = 0)
    : BranchPlain(b)
    , _height{h}
  {}
  virtual ~NodeAge(){}

  void setHeight(double height) { _height = height; }
  double getHeight(void) const { return _height; }

  // TODO overwirte serialize/deserialize

private: 
  double _height;  
};




#endif /* NODEAGE_H */
