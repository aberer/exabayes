#ifndef _BLOCK_PARTITION_H
#define _BLOCK_PARTITION_H

#include <ncl/ncl.h>

#include "GlobalVariables.hpp"
#include "RandomVariable.hpp"
#include "TreeAln.hpp"

class BlockParams : public NxsBlock
{
public: 
  explicit BlockParams(int numPart)
    : numPart(numPart)
    , hasAA(false)
    , hasDna(false)
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PARAMS";    

    
    // TODO use hasAA! 
  }

  vector<RandomVariable> getParameters() const{return parameters; }
  virtual void Read(NxsToken &token); 
  void initialize(const TreeAln &traln); 


private: 
  void parseScheme(NxsToken& token, category_t cat, nat &idCtr); 
  

  int numPart;   
  vector<RandomVariable> parameters; 

  bool hasAA;		 // does our alignment have an aa partition? 
  bool hasDna;			// does our alignment have a dna partition?  
}; 


#endif
