
#ifndef _UTILS_TOPO_H
#define _UTILS_TOPO_H



void insertWithGenericBL (nodeptr insertNode, nodeptr branchNode, double *insertZ, double *branchNodeZ, double *neighbourZ ,  int numBranches);
void insertWithUnifBL (nodeptr insertNode, nodeptr branchNode, int numBranches);
void insertWithUnifBLScaled(nodeptr insertNode, nodeptr branchNode, double scale, int numBranches);

void traverseAndCount(nodeptr p, int *count, tree *tr ); 
#endif
