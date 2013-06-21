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
  virtual void applyToTree(TreeAln &traln) const ; 
  virtual void revertTree(TreeAln &traln, PriorBelief &prior) const ; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const ; 
  virtual void extractMoveInfo(const TreeAln &traln, vector<Branch> description) ; 
  virtual void multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, vector<shared_ptr<AbstractPrior> > brPr) const ; 

  virtual Branch getEvalBranch(const TreeAln &traln) const {return innerBranch; }
  virtual AbstractMove* clone() const {return new NniMove();}  

private: 			// METHODS 
  void multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter , vector<shared_ptr<AbstractPrior> > priors, string name) const; 
  
private: 			// ATTRIBUTES 
  Branch innerBranch; 
  vector<Branch> outerBranches; 
  pair<int,int> switching; 
}; 



#endif
