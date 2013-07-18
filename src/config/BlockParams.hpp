

#ifndef _BLOCK_PARTITION_H
#define _BLOCK_PARTITION_H

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "parameters/AbstractParameter.hpp" 
#include "TreeAln.hpp"
#include "Category.hpp"

class BlockParams : public NxsBlock
{
public: 
  BlockParams()
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PARAMS";    
  }

  void setTree(shared_ptr<TreeAln> _traln){traln = _traln; }
  vector<unique_ptr<AbstractParameter> > getParameters() const; 
  virtual void Read(NxsToken &token); 

private: 
  void parseScheme(NxsToken& token, Category cat, nat &idCtr); 
  vector<unique_ptr<AbstractParameter> > parameters; 
  shared_ptr<TreeAln> traln; 
}; 


#endif
