#ifndef _BLOCK_PRIOR_H
#define _BLOCK_PRIOR_H

#include <map>
#include <memory>
#include <unordered_map>

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "Priors.hpp"
#include "axml.h"

#include "Category.hpp"

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

  vector< PriorPtr > getGeneralPriors() const {return generalPriors; }
  vector< map<nat,PriorPtr>>  getSpecificPriors() const {return specificPriors; }

private: 
  nat numPart;   
  vector<PriorPtr> generalPriors; 
  vector< map<nat,PriorPtr>> specificPriors; // for each category
}; 


#endif
