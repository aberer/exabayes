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




void TreePrinter::helper(const TreeAln &traln, stringstream &ss, nodeptr p)
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
      ss << "("; 
      helper(traln, ss, p->next->back); 
      ss << "," ; 
      helper(traln, ss, p->next->next->back);       
      ss << ")"; 
    }
  
  if(p == traln.getTr()->start->back)    
    {
      ss << ";"; 
      return; 
    }
  
  if(withInternalNodes && not traln.isTipNode(p))
    ss << p->number << "[x:" << getTheX(p) << "]"; 

  if(withBranchLengths)
    ss << ":" << setprecision(7) << fixed  << branchLengthToReal(traln.getTr(), p->z[0]); 
}

string TreePrinter::printTree(const TreeAln &traln)
{
  stringstream ss; 
  helper(traln, ss, traln.getTr()->start->back);   
  return ss.str(); 
}
