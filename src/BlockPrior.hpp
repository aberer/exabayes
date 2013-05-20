#ifndef _BLOCK_PRIOR_H
#define _BLOCK_PRIOR_H

#include <map>
#include <memory>
#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "Priors.hpp"
#include "axml.h"


class BlockPrior : public NxsBlock
{
public: 
  explicit BlockPrior(nat numPart) 
    : numPart(numPart) 
    , generalPriors(NUM_PROP_CATS)
    , specificPriors(NUM_PROP_CATS)
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PRIOR"; 
  }
 

  shared_ptr<AbstractPrior> parsePrior(NxsToken &token)  ; 
  virtual void Read(NxsToken &token); 

private: 
  nat numPart; 
  
  vector<shared_ptr<AbstractPrior>> generalPriors; 

  vector< map<nat, shared_ptr<AbstractPrior> > > specificPriors; // for each category
}; 


#endif
