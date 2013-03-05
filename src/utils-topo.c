#include "globals.h"
#include "common.h"

#include "config.h"
#include "axml.h"

#include "main-common.h"


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
