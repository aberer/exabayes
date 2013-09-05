#ifndef _TBR_MOVE_H 
#define _TBR_MOVE_H 

#include "SprMove.hpp"

class TbrMove : public SprMove
{
public: 
  virtual ~TbrMove(){}   
  virtual void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const ; 
  virtual void revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const ; 
  virtual void disorientAtNode(TreeAln &traln, nodeptr p) const ; 
  virtual void extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description, const std::vector<AbstractParameter*> &params) ; 
  virtual BranchPlain getEvalBranch(const TreeAln &traln) const ; 
  virtual AbstractMove* clone() const {return new TbrMove;}  

  friend std::ostream& operator<<(std::ostream &out, const TbrMove &rhs) ;  

private: 			// ATTRIBUTES
  Path path1; 
  Path path2; 
}; 

#endif
