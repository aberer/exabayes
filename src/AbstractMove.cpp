#include "AbstractMove.hpp"

void AbstractMove::disorientHelper(const TreeAln &traln, nodeptr p) const 
{
  if(traln.isTipNode(p))
    {
    }
  else if(p->x)
    {
      p->x = 0; 
      p->next->x = 1; 
    }
  else 
    {
    }
}
