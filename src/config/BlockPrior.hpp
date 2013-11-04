#ifndef _BLOCK_PRIOR_H
#define _BLOCK_PRIOR_H

#include <memory>
#include <unordered_map>
#include <unordered_set>

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "priors/AbstractPrior.hpp"
#include "axml.h"

#include "Category.hpp"


// if the set is empty, then we have a general "fall-back" prior
typedef std::unordered_multimap<Category, 
				std::tuple<std::unordered_set<nat>,
					   std::unique_ptr<AbstractPrior> > >  
multiMapCategory2TuplePartitionsPrior ; 


class BlockPrior : public NxsBlock
{
public: 
  explicit BlockPrior(nat numPart) 
    : _numPart(numPart)
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PRIOR"; 
  }
  
  virtual void Read(NxsToken &token); 
  const multiMapCategory2TuplePartitionsPrior& getPriors()const  {return _parsedPriors; } 

private: 			// METHODS
  std::unique_ptr<AbstractPrior> parsePrior(NxsToken &token)  ; 
  
private: 			// ATTRIBUTES
  multiMapCategory2TuplePartitionsPrior _parsedPriors; 
  nat _numPart;
}; 


#endif
