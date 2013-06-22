

#ifndef _BLOCK_PARTITION_H
#define _BLOCK_PARTITION_H

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "RandomVariable.hpp"
#include "TreeAln.hpp"

class BlockParams : public NxsBlock
{
public: 
  explicit BlockParams(const TreeAln& _traln)
    :traln(_traln)
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PARAMS";    
  }

  vector<RandomVariable> getParameters() const{return parameters; }
  virtual void Read(NxsToken &token); 
  void initialize(const TreeAln &traln); 

private: 
  void parseScheme(NxsToken& token, category_t cat, nat &idCtr); 
  vector<RandomVariable> parameters; 
  const TreeAln &traln; 
}; 


#endif
