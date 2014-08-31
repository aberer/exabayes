#ifndef NODEAGE_H
#define NODEAGE_H

// for this class, _primNode is the node, we identify the age for 

// _secNode points upwards (with the NOT having any upward node)


class NodeAge : public BranchPlain
{
public:
  NodeAge(){}
  virtual ~NodeAge(){}
  
  
private: 
  double _height;  
};




#endif /* NODEAGE_H */
