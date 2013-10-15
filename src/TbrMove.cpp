#include "Branch.hpp"
#include "TbrMove.hpp"
#include "Path.hpp"

void TbrMove::applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const 
{
  applyPath(traln,path1, blParams); 
  applyPath(traln,path2, blParams); 
} 


void TbrMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const 
{
  auto path1r = getPathAfterMove(path1);
  auto path2r = getPathAfterMove(path2);
  
  applyPath(traln,path1r, blParams);
  applyPath(traln,path2r, blParams);

  path1.restoreBranchLengthsPath(traln, blParams); 
  path2.restoreBranchLengthsPath(traln, blParams); 
} 


void TbrMove::extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description,const std::vector<AbstractParameter*> &params) 
{
  sprCreatePath(traln, description.at(0),description.at(1),path1, params); 
  sprCreatePath(traln, description.at(2),description.at(3),path2, params); 
} 


BranchPlain TbrMove::getEvalBranch(const TreeAln &traln) const 
{   
  return  BranchPlain(path1.getNthNodeInPath(1) , path2.getNthNodeInPath(1)); 
}


std::ostream& operator<<(std::ostream &out, const TbrMove &rhs) 
{
  return out <<  "path1:"  << rhs.path1 << std::endl
	     << "path2:" << rhs.path2 << std::endl; 
}  


std::vector<nat> TbrMove::getDirtyNodes(const TreeAln &traln, bool considerOuter) const
{
  assert(not considerOuter); 

  auto result = std::vector<nat>(); 
  result.reserve(path1.getNumberOfNodes() + path2.getNumberOfNodes() );
  
  for(int i = 1; i < path1.getNumberOfNodes()- 1; ++i)
    result.push_back ( path1.getNthNodeInPath(i)) ;

  for(int i = 1 ; i < path2.getNumberOfNodes() - 1 ; ++i)
    result.push_back( path2.getNthNodeInPath(i)) ; 

  return result; 
} 

