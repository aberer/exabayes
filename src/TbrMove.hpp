#ifndef _TBR_MOVE_H 
#define _TBR_MOVE_H 

#include "SprMove.hpp"

class TbrMove : public SprMove
{
public: 
  void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const ; 
  void revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const ; 
  void extractMoveInfo(const TreeAln &traln, std::tuple<BranchPlain,BranchPlain,BranchPlain> description,const std::vector<AbstractParameter*> &params) ; 
  BranchPlain getEvalBranch(const TreeAln &traln) const ; 
  std::vector<nat> getDirtyNodes(const TreeAln &traln, bool considerOuter) const; 

  friend std::ostream& operator<<(std::ostream &out, const TbrMove &rhs) ;  

private: 			// ATTRIBUTES
  Path path1; 
  Path path2; 
}; 

#endif
