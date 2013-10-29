#ifndef _BLOCK_PRIOR_H
#define _BLOCK_PRIOR_H

#include <memory>
#include <unordered_map>

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "priors/AbstractPrior.hpp"
#include "axml.h"

#include "Category.hpp"

class BlockPrior : public NxsBlock
{
public: 
  explicit BlockPrior(nat numPart) 
    : numPart(numPart) 
    , generalPriors(CategoryFuns::getAllCategories().size())
    , specificPriors(CategoryFuns::getAllCategories().size())
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PRIOR"; 
  }
 

  shared_ptr<AbstractPrior> parsePrior(NxsToken &token)  ; 
  virtual void Read(NxsToken &token); 

  vector< shared_ptr<AbstractPrior> > getGeneralPriors() const {return generalPriors; }
  std::vector< std::unordered_map<nat,std::shared_ptr<AbstractPrior>>>  getSpecificPriors() const {return specificPriors; }

private: 
  nat numPart;   
  std::vector<shared_ptr<AbstractPrior>> generalPriors; 
  std::vector< std::unordered_map<nat,std::shared_ptr<AbstractPrior>>> specificPriors; // for each category
}; 


#endif
