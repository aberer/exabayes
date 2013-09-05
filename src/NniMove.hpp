#ifndef _NNI_MOVE
#define _NNI_MOVE

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "PriorBelief.hpp"

#include "AbstractMove.hpp"

class  NniMove : public AbstractMove
{
public: 
  virtual void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const ; 
  virtual void revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const ; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const ; 
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description, const std::vector<AbstractParameter*> &params) ; 
  // virtual void multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<AbstractPrior* > brPr) const ; 
  virtual BranchPlain getEvalBranch(const TreeAln &traln) const { return innerBranch.toPlain(); }
  virtual AbstractMove* clone() const {return new NniMove();}  

  friend std::ostream& operator<<(std::ostream& out, const NniMove &move); 

  /**
     @brief maps a branch that exists in the tree after executing the
     move to the respective branch in the tree before the move was
     executed.
   */ 
  template<typename BT>
  BranchPlain mapToBranchBeforeMove(const BT &b) const ; 
  
private: 			// METHODS 
  // void multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter , std::vector<AbstractPrior* > priors, std::string name) const; 
  
private: 			// ATTRIBUTES 
  BranchLengths innerBranch; 
  std::vector<BranchLengths> outerBranches; 
  std::pair<int,int> switching; 
}; 



#endif
