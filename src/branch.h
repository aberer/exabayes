#ifndef _BRANCH_H
#define _BRANCH_H


struct _state;  


typedef struct 
{
  int thisNode;  /// in case we want to give the branch some orientation, this is the node we want to extract   
  int thatNode;  

  /* todo also branch lengths?  */
} branch; /// represents a branch 


void insertNodeIntoBranch(struct _state  *chain, branch toBeInserted, branch insertionBranch, double* blsNode1, double* blsNode2); 
branch invertBranch(branch b); 
branch constructBranch(int thisNode, int thatNode); 
boolean branchExists(tree *tr, branch b); 
nodeptr findNodeFromBranch(tree *tr, branch b ); 
void pruneBranch(struct _state  *chain, branch b, double *z); 
branch findRoot(struct _state *chain); 
#endif
