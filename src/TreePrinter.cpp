#include <cassert>
#include <iomanip>

#include "Branch.hpp"
#include "TreePrinter.hpp"

// TODO replace with ostream 


TreePrinter::TreePrinter(bool withBranchLengths , bool withInternalNodes , bool withRealNames) 
  : withBranchLengths(withBranchLengths)
  , withInternalNodes(withInternalNodes)
  , withRealNames(withRealNames)
{
}


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


std::string TreePrinter::printTree(const TreeAln& traln, const std::vector<AbstractParameter*> &params)
{
  std::stringstream ss; 
  ss << SOME_FIXED_PRECISION; 

  ss << "("; 
  helper(traln, ss, traln.getTr()->start->back, true, params );   
  ss << ")" ; 
  if(withBranchLengths)
    ss << ":0.0"; 
  ss << ";"; 
  return ss.str();
} 


void TreePrinter::printBranchLength(const TreeAln& traln, std::stringstream &ss, nodeptr p , const std::vector<AbstractParameter*> &params)
{
  ss << MAX_SCI_PRECISION; 

  if(params.size() == 1 )
    {
      auto param = params[0]; 
      ss << ":" << traln.getBranch(p,param).getInterpretedLength(traln, param); 
    }
  else 
    {
      ss << ":["; 
      bool isFirst = true ; 
      for(auto &param : params)
	{
	  ss << (isFirst ? "" : ",")  << traln.getBranch(p,param).getInterpretedLength(traln, param); 
	  isFirst = false; 
	}
      ss << "]"; 
    }
}


void TreePrinter::helper(const TreeAln &traln, std::stringstream &ss, 
			 nodeptr p, bool isFirst, const std::vector<AbstractParameter*> &params)
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
      helper(traln, ss, p->next->back, false, params); 
      ss << "," ; 
      helper(traln, ss, p->next->next->back, false, params ); 
      if(not isFirst)
	ss << ")"; 
    }

  if(not isFirst && withInternalNodes && not traln.isTipNode(p))
    ss << p->number ; 

  if(not isFirst && withBranchLengths && params.size() != 0 )
    printBranchLength(traln, ss, p, params);

  if(isFirst)
    {
      ss << "," << p->back->number; 
      if(withBranchLengths && params.size() !=  0)
	printBranchLength(traln, ss, p, params);
    }
}

std::string TreePrinter::printTree(const TreeAln& traln)
{
  return printTree(traln, {});
} 
