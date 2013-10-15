#ifndef _ABSTRACT_MOVE_H
#define _ABSTRACT_MOVE_H

#include <vector>
#include <memory>
#include "BranchFwd.hpp"
#include "axml.h"

class TreeAln; 
class PriorBelief; 
class AbstractPrior; 
class Path; 
class AbstractParameter; 

class AbstractMove
{
public: 
  virtual AbstractMove* clone() const = 0; 
  virtual ~AbstractMove(){}
  
  virtual void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const = 0; 
  virtual void revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const = 0; 
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description, const std::vector<AbstractParameter*> &params ) = 0; 
  virtual BranchPlain getEvalBranch(const TreeAln &traln) const = 0; 


  virtual std::vector<nat> getDirtyNodes(const TreeAln &traln, bool considerOuter) const = 0; 
}; 


#endif
