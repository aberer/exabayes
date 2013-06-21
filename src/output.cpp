#include <sstream>
#include <cassert>

#include "axml.h"
#include "GlobalVariables.hpp"
#include "output.h"
#include "TreeAln.hpp"
#include "branch.h"



/**
   @brief A check to ensure that the tree is still consistent. 
   
   @param  p -- pointer to an inner (!) node
 */ 
static void traverseAndCount(nodeptr p, int &count, const tree *tr )
{
  nodeptr q;  

  ++count;

  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseAndCount(q->back, count, tr);
	  q = q->next;
	} 
    }
}



/**
   @brief prints the tree as a debug message. 
*/
void debug_printTree(TreeAln &traln)
{
#ifdef DEBUG_SHOW_TREE
  TreePrinter tp(true, false, false); 
  tout << "tree: " << tp.printTree(traln) << endl; 
#endif
}


void debug_checkTreeConsistency(const TreeAln &traln)
{
#ifdef DEBUG_CHECK_TREE_CONSISTENCY
  auto tr = traln.getTr(); 
  int count = 0; 
  traverseAndCount(tr->start->back, count, tr); 
  if(count != 2 * tr->mxtips - 3 )
    { 
      TreePrinter tp(true, false, false); 
      string treeString = tp.printTree(traln); 

      tout << "faulty Topology: " << treeString << endl; 

      assert(2 * tr->mxtips-3 == count); 
    }
#endif
}



void printOrientation(tree *tr, nodeptr p)
{
  if(isTip(p->number,tr->mxtips))
    printf("%d is tip\n", p->number); 
  else if(p->x )
    printf("%d is orientated towards %d\n", p->number,p->back->number); 
  else if(p->next->x)
    printf("%d is orientated towards %d\n", p->number,p->next->back->number); 
  else if(p->next->next->x)
    printf("%d is orientated towards %d\n", p->number,p->next->next->back->number); 
  else 
    assert(0); 
}

