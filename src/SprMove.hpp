/** 
    @brief a helper object making the code concerning a SPR move
    re-usable.
    
    @todo still could use some cleanup.
 */ 

#ifndef _SPR_MOVE_H
#define _SPR_MOVE_H

#include "axml.h"
#include "AbstractMove.hpp"
#include "Path.hpp"

class PriorBelief; 
class LikelihoodEvaluator; 


class SprMove : public AbstractMove
{
public:    		
  virtual void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const; 
  virtual void revertTree(TreeAln &traln,const std::vector<AbstractParameter*> &blParams ) const; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const; 
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description,const std::vector<AbstractParameter*> &params); 
  virtual BranchPlain getEvalBranch(const TreeAln &traln) const; 

  std::vector<BranchLength> proposeBranches(TreeAln &traln, const std::vector<AbstractParameter*> &blParams, LikelihoodEvaluator& eval, double &hastings, Randomness& rand, bool isForward); 
  virtual AbstractMove* clone() const; 
  /** 
      @brief gets the number of nni moves to be executed in order to
      achieve this move
  */ 
  nat getNniDistance() const { return path.size( )- 2 ;  }

  void integrateBranches( TreeAln &traln,  const std::vector<AbstractParameter*> blParam, LikelihoodEvaluator &eval, double &hastings ) const ; 

  // std::vector<BranchLength> proposeBranches(TreeAln &traln); 

  BranchPlain getInsertionBranchAfterInner() const ; 
  BranchPlain getInsertionBranchAfterOuter()  const ; 
  BranchPlain getInsertionBranchBefore() const ; 
  BranchPlain getSubtreeBranchBefore(const TreeAln &traln ) const ; 
  BranchPlain getSubtreeBranchAfter(const TreeAln &traln ) const ; 
  BranchPlain getPruningBranchBeforeInner() const ; 
  BranchPlain getPruningBranchBeforeOuter() const; 
  BranchPlain getPruningBranchAfterPrune()const ;   
  BranchPlain getOppositeBranch(const TreeAln &traln ) const ; 

  friend std::ostream& operator<<(std::ostream &out, const SprMove& rhs); 

protected:			// METHODS
  void sprCreatePath(const TreeAln &traln, BranchPlain mover, BranchPlain movedInto, Path &path, const std::vector<AbstractParameter*> &params ) const;
  void sprDisorientPath(TreeAln &traln, nodeptr p, const Path &path) const ;
  void applyPath(TreeAln &traln, const Path &modifiedPath, const std::vector<AbstractParameter*> &params ) const; 
  void getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath) const; 

private: 			// ATTRIBUTES
  Path path; 
}; 

#endif
