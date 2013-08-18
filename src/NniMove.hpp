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
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<Branch> description, const std::vector<AbstractParameter*> &params) ; 
  // virtual void multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<AbstractPrior* > brPr) const ; 
  virtual Branch getEvalBranch(const TreeAln &traln) const {return innerBranch; }
  virtual AbstractMove* clone() const {return new NniMove();}  

  virtual NniMove getInvertseMove() const ; 

  
private: 			// METHODS 
  // void multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter , std::vector<AbstractPrior* > priors, std::string name) const; 
  
private: 			// ATTRIBUTES 
  Branch innerBranch; 
  std::vector<Branch> outerBranches; 
  std::pair<int,int> switching; 
}; 



#endif
