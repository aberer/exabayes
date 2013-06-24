

#ifndef _BLOCK_PARTITION_H
#define _BLOCK_PARTITION_H

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "RandomVariable.hpp"
#include "TreeAln.hpp"
#include "Category.hpp"

class BlockParams : public NxsBlock
{
public: 
  BlockParams()
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PARAMS";    
  }

  void setTree(TreeAlnPtr _traln){traln = _traln; }
  vector<RandomVariablePtr> getParameters() const{return parameters; }
  virtual void Read(NxsToken &token); 

private: 
  void parseScheme(NxsToken& token, Category cat, nat &idCtr); 
  vector<RandomVariablePtr> parameters; 
  TreeAlnPtr traln; 
}; 


#endif
