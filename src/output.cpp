#include <sstream>
#include <cassert>

#include "axml.h"
#include "GlobalVariables.hpp"
#include "output.h"		// 
// #include "adapters.h"
#include "TreeAln.hpp"
#include "branch.h"

bool isOutputProcess()
{
#if HAVE_PLL != 0
  return true; 
#else 
  return processID == 0; 
#endif   
}


/**
   @brief A check to ensure that the tree is still consistent. 
   
   @param  p -- pointer to an inner (!) node
 */ 
static void traverseAndCount(nodeptr p, int *count, tree *tr )
{
  nodeptr q;  

  *count += 1;

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
  tree *tr = traln.getTr();
  Tree2stringNexus(tr->tree_string, tr, tr->start->back, 0); 
  tout << "tree: "<< tr->tree_string << endl; 
#endif
}




/**
   @brief prints the entire environment of a node with branch lengths
 */
void debug_printNodeEnvironment(TreeAln &traln, int nodeID )
{
#ifdef DEBUG_SHOW_TREE
  tree *tr = traln.getTr();
  nodeptr 
    p =  tr->nodep[nodeID]; 

  
  printf("STATE node %d:\thooked to %d (bl=%f),\t%d (bl=%f)\tand %d (bl=%f)\n", p->number, p->back->number, traln.getBranchLength( p->back,0), 
	 p->next->back->number, traln.getBranchLength( p->next->back,0),
	    p->next->next->back->number, traln.getBranchLength( p->next->next->back,0)
	    );   
#endif
}


void debug_checkTreeConsistency(tree *tr)
{
#ifdef DEBUG_CHECK_TREE_CONSISTENCY
  // tree *tr = chain->tr; 
  int count = 0; 
  traverseAndCount(tr->start->back, &count, tr); 
  if(count != 2 * tr->mxtips - 3 )
    {      
      char tmp[10000]; 
      Tree2stringNexus(tmp, tr, tr->start->back, 0); 
      if(isOutputProcess())
	printf("faulty TOPOLOGY: %s\n", tmp);

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


char *Tree2stringNexus(char *treestr, tree *tr , nodeptr p, int perGene )
{   
  // tree *tr = chain->tr; 

  if(isTip(p->number, tr->mxtips)) 
    {	       	        
      sprintf(treestr, "%d", p->number); 
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = Tree2stringNexus(treestr, tr, p->next->back, perGene);
      *treestr++ = ',';
      treestr = Tree2stringNexus(treestr, tr, p->next->next->back, perGene);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2stringNexus(treestr, tr, p->back, perGene);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back) 
    {
      sprintf(treestr, ";"); 
      return treestr;
    }

  sprintf(treestr, ":%.7f", branchLengthToReal(tr, p->z[0]));
  
  while (*treestr) treestr++;
  return  treestr;
}
