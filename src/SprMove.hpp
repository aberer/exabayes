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
  virtual void multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, vector<shared_ptr<AbstractPrior> > brPr ) const; 
  virtual Branch getEvalBranch(const TreeAln &traln) const; 

  virtual AbstractMove* clone() const; 

  friend ostream& operator<<(ostream &out, const SprMove& rhs); 

protected:			// METHODS
  void sprCreatePath(const TreeAln &traln, Branch mover, Branch movedInto, Path &path ) const;
  void sprDisorientPath(TreeAln &traln, nodeptr p, const Path &path) const ;
  void applyPath(TreeAln &traln, const Path &modifiedPath ) const; 
  void getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath) const; 

private: 			// ATTRIBUTES
  Path path; 
}; 

#endif
