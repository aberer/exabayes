

#ifndef _BLOCK_PARTITION_H
#define _BLOCK_PARTITION_H

#include "config/ExaBlock.hpp"

#include "system/GlobalVariables.hpp"
#include "parameters/AbstractParameter.hpp" 
#include "model/TreeAln.hpp"
#include "model/Category.hpp"

class BlockParams : public ExaBlock
{
public: 
  BlockParams()
  {
    NCL_BLOCKTYPE_ATTR_NAME = "PARAMS";    
  }

  void setTree(const TreeAln* _traln){ tralnPtr = _traln; }
  vector<unique_ptr<AbstractParameter> > getParameters() const; 
  virtual void Read(NxsToken &token); 

private:   			// METHODS
  void partitionError(nat partition, nat totalPart) const ; 
  void parseScheme(NxsToken& token, Category cat, nat &idCtr); 
  nat getNumSeen(Category cat) ; 
  
private: 			// ATTRIBUTES
  vector<unique_ptr<AbstractParameter> > parameters; 
  const TreeAln* tralnPtr;  	// NON-owning
}; 


#endif
