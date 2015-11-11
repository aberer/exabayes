#include "model/Branch.hpp"
#include "TbrMove.hpp"
#include "data-struct/Path.hpp"


void TbrMove::applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const 
{
  if( path1.size() > 0 )
    applyPath(traln,path1, blParams); 
  if( path2.size() > 0)
    applyPath(traln,path2, blParams); 
} 


void TbrMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const 
{
  if( path1.size() > 0)
    {
      auto path1r = getPathAfterMove(path1);
      applyPath(traln,path1r, blParams);
      path1.restoreBranchLengthsPath(traln, blParams); 
    }

  if( path2.size() > 0)
    {
      auto path2r = getPathAfterMove(path2);
      applyPath(traln,path2r, blParams);
      path2.restoreBranchLengthsPath(traln, blParams); 
    }
} 


void TbrMove::extractMoveInfo(const TreeAln &traln, std::tuple<BranchPlain,BranchPlain,BranchPlain> description,const std::vector<AbstractParameter*> &params) 
{
  path1.clear(); 
  path2.clear(); 
  
  auto 
    empty = BranchPlain{},
    elem0 = std::get<0>(description),
    elem1 = std::get<1>(description),
    elem2 = std::get<2>(description); 

  if(not elem1.equalsUndirected(empty))
    sprCreatePath(traln, elem0,elem1,path1, params); 
  if(not elem2.equalsUndirected(empty))
    sprCreatePath(traln, elem0.getInverted(),elem2,path2, params); 
} 


BranchPlain TbrMove::getEvalBranch(const TreeAln &traln) const 
{   
  auto evalBranch = BranchPlain{}; 

  if(path1.size( ) > 0 && path2.size( ) ) 
    evalBranch = BranchPlain(path1.getNthNodeInPath(1) , path2.getNthNodeInPath(1)); 
  else if(path1.size( ) > 0 )
    evalBranch = getEvalBranchFromPath(traln, path1); 
  else if(path2.size() > 0 )
    evalBranch = getEvalBranchFromPath(traln, path2); 
  else 
    assert(0); 

  return evalBranch; 
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

  if(path1.size( ) > 0 )
    {
      for(int i = 1; i < path1.getNumberOfNodes()- 1; ++i)
	result.push_back ( path1.getNthNodeInPath(i)) ;
    }

  if(path2.size() > 0 )
    {
      for(int i = 1 ; i < path2.getNumberOfNodes() - 1 ; ++i)
	result.push_back( path2.getNthNodeInPath(i)) ; 
    }

  return result; 
} 

