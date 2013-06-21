/** 
    @brief a helper object making the code concerning a SPR move
    re-usable.
    
    @todo still could use some cleanup.
 */ 

#ifndef _SPR_MOVE_H
#define _SPR_MOVE_H

#include "axml.h"
#include "PriorBelief.hpp"
#include "Path.hpp"
#include "AbstractMove.hpp"


class SprMove : public AbstractMove
{
public:    
  virtual void applyToTree(TreeAln &traln) const; 
  virtual void revertTree(TreeAln &traln, PriorBelief &prior) const; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const; 
  virtual void extractMoveInfo(const TreeAln &traln, vector<Branch> description); 
  virtual void multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, shared_ptr<AbstractPrior> brPr) const; 
  virtual Branch getEvalBranch(const TreeAln &traln) const; 

  virtual AbstractMove* clone() const; 

  friend ostream& operator<<(ostream &out, const SprMove& rhs); 

  // legacy 
  void multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior,  const Path &modifiedPath, double multiplier, shared_ptr<AbstractPrior> brPr) const; 
  void applyPathAsESPR(TreeAln &traln, const Path &modifiedPath ) const; 
  void resetAlongPathForESPR(TreeAln &traln, PriorBelief &prior , const Path &modifiedPath) const;   
  void getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath) const; 




private: 
  Path path; 
}; 

#endif
