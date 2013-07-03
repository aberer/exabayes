#include <vector>

#include "TreeRandomizer.hpp"
#include "TreeAln.hpp"
#include "treeRead.h"


TreeRandomizer::TreeRandomizer(randCtr_t seed)
  : rand(seed)
{
}

void TreeRandomizer::randomizeTree(TreeAln &traln)
{
  for(nat i = 1 ; i < traln.getNumberOfNodes() + 1 ; ++i)
    {
      auto p = traln.getNode(i);
      p->back = NULL; 
      p->next->back = NULL; 
      p->next->next->back = NULL;       
    }

  // start with the simple tree 
  auto a = traln.getNode(1 ),
    b = traln.getNode(2 ), 
    c = traln.getNode( 3 ), 
    inner = traln.getNode(traln.getNumberOfTaxa() + 1); 
  
  traln.clipNodeDefault(inner,a); 
  traln.clipNodeDefault(inner->next, b); 
  traln.clipNodeDefault(inner->next->next,c); 

  for(nat i = 4; i < traln.getNumberOfTaxa() +1 ;  ++i)
    {
      inner = traln.getNode(traln.getNumberOfTaxa() + i-2 );       
      auto taxonP = traln.getNode(i); 
      traln.clipNodeDefault(taxonP, inner); 
      
      auto b = traln.drawBranchUniform_helper(rand, i-1);
      auto p1 = b.findNodePtr(traln),
	p2 = p1->back; 
      
      traln.clipNodeDefault(p1, inner->next); 
      traln.clipNodeDefault(p2, inner->next->next); 
    }
  
  traln.getTr()->start = traln.getTr()->nodep[1];   
}
