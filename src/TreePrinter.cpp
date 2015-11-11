#include "TreePrinter.hpp"

static int getTheX(nodeptr p)
{
  if(p->x)
    return p->back->number; 
  else if(p->next->x)
    return p->next->back->number; 
  else if(p->next->next->x)
    return p->next->next->back->number; 
  else 
    assert(0) ; 
  return 0;     
}




void TreePrinter::helper(const TreeAln &traln, std::stringstream &ss, nodeptr p, bool isFirst)
{
  if(traln.isTipNode(p))
    {
      if(withRealNames)
	ss << traln.getTr()->nameList[p->number]; 
      else 
	ss << p->number; 
    }
  else 
    {
      if (not isFirst )
	ss << "("; 
      helper(traln, ss, p->next->back, false); 
      ss << "," ; 
      helper(traln, ss, p->next->next->back, false );       
      if(not isFirst)
	ss << ")"; 
    }

  if(not isFirst && withInternalNodes && not traln.isTipNode(p))
    ss << p->number ; 

  if(not isFirst && withBranchLengths)
    ss << ":" << std::setprecision(7) << std::fixed  << Branch(p->number,p->back->number,p->z[0]).getInterpretedLength(traln); 

  if(isFirst)
    {
      ss << "," << p->back->number; 
      if(withBranchLengths)
	ss << ":" << std::setprecision(7) << std::fixed  << Branch(p->back->number, p->number, p->z[0]).getInterpretedLength(traln); 
    }
}


// this is unrooted! 
std::string TreePrinter::printTree(const TreeAln &traln)
{
  std::stringstream ss; 
  ss << "("; 
  helper(traln, ss, traln.getTr()->start->back, true );   
  ss << "):0.0;"; 

  return ss.str();

}
