/**
   @file  topology-utils.c
   
   @brief Various functions for manipulation of topologies and branch lengths. 
*/ 
#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "randomness.h"
#include "main-common.h"
#include "adapters.h"
#include "output.h"



/**
   @brief A check to ensure that the tree is still consistent. 
   
   @param  p -- pointer to an inner (!) node
 */ 
void traverseAndCount(nodeptr p, int *count, tree *tr )
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




void insertWithGenericBL (nodeptr insertNode, nodeptr branchNode, double *insertZ, double *branchNodeZ, double *neighbourZ ,  int numBranches)
{
  nodeptr thirdNode = branchNode->back;   

  /* if(processID == 0) */
  /*   { */
  /*     printf("insertNode=%d, branchNode=%d, thirdNode=%d\n", insertNode->number, branchNode->number, thirdNode->number);  */
  /*   } */

  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
  hookup(insertNode, insertNode->back, insertZ, numBranches);
}

void insertWithUnifBL (state* chain, nodeptr insertNode, nodeptr branchNode, int numBranches)
{
  /* BUG: drawRandDouble does not take an argument  */
  assert(0); 

  double r;
  nodeptr thirdNode = branchNode->back;
  double branchNodeZ[numBranches];
  double neighbourZ[numBranches];
  
  for(int i=0; i<numBranches;i++)
  {
   /* r=drawRandDouble(branchNode->z[i]); */
   branchNodeZ[i]=r;
   neighbourZ[i]=1-r; 
  }
  
  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
}



void insertWithUnifBLScaled(nodeptr insertNode, nodeptr branchNode, double scale, int numBranches)
{
  /* BUG: drawRandDouble does not take an argument  */
  assert(0); 


  double r;
  nodeptr thirdNode = branchNode->back;
  double branchNodeZ[numBranches];
  double neighbourZ[numBranches];
  
  for(int i=0; i<numBranches;i++)
  {
   /* r=drawRandDouble(scale*branchNode->z[i]); */
   branchNodeZ[i]=r;
   neighbourZ[i]=1-r; 
  }
  
  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
}

