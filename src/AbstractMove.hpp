#ifndef _ABSTRACT_MOVE_H
#define _ABSTRACT_MOVE_H

#include <memory>

#include "TreeAln.hpp"
#include "PriorBelief.hpp"
#include "Branch.hpp"
#include "priors/AbstractPrior.hpp"
#include "Path.hpp"


class AbstractMove
{
public: 
  virtual AbstractMove* clone() const = 0; 
  virtual ~AbstractMove(){}
  
  virtual void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const = 0; 
  virtual void revertTree(TreeAln &traln, PriorBelief &prior, const std::vector<AbstractParameter*> &params) const = 0; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const = 0; 
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<Branch> description, const std::vector<AbstractParameter*> &params ) = 0; 
  virtual Branch getEvalBranch(const TreeAln &traln) const = 0; 

protected: 
  void disorientHelper(const TreeAln &traln, nodeptr p) const ; 


}; 


#endif
