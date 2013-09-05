#include "AbstractMove.hpp"
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"

void AbstractMove::disorientHelper(const TreeAln &traln, nodeptr p) const 
{
  if(traln.isTipNode(p))
    {
    }
  else if(p->x)
    {
      // tout << "disorienting " << p->number << std::endl; 
      p->x = 0; 
      p->next->x = 1; 
    }
  else 
    {
      // tout << "already " << p->number << std::endl; 
    }
}

