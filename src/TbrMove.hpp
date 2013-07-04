#ifndef _TBR_MOVE_H 
#define _TBR_MOVE_H 

#include "SprMove.hpp"

class TbrMove : public SprMove
{
public: 
  virtual ~TbrMove(){}   
  virtual void applyToTree(TreeAln &traln) const ; 
  virtual void revertTree(TreeAln &traln, PriorBelief &prior) const ; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const ; 
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<Branch> description) ; 
  virtual void multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<AbstractPrior*> brPr) const ; 
  virtual Branch getEvalBranch(const TreeAln &traln) const ; 
  virtual AbstractMove* clone() const {return new TbrMove;}  

  friend std::ostream& operator<<(std::ostream &out, const TbrMove &rhs) ;  

private: 			// ATTRIBUTES
  Path path1; 
  Path path2; 
}; 

#endif
