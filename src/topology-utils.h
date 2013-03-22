
#ifndef _UTILS_TOPO_H
#define _UTILS_TOPO_H

nodeptr findNodePtrByNeighbor(tree *tr, int targetNode,  int neighborNode); 

void insertWithGenericBL (nodeptr insertNode, nodeptr branchNode, double *insertZ, double *branchNodeZ, double *neighbourZ ,  int numBranches);
void insertWithUnifBL (nodeptr insertNode, nodeptr branchNode, int numBranches);
void insertWithUnifBLScaled(nodeptr insertNode, nodeptr branchNode, double scale, int numBranches);


void pruneNodeFromNodes(state *chain , int prunedNode,  int neighbor, int nextNeighbor, double *bls); 
void insertNodeIntoBranch(state *chain, int insertNode, int branchNode1, int branchNode2, double* blsNode1, double* blsNode2); 


void traverseAndCount(nodeptr p, int *count, tree *tr ); 

boolean nodesAreHooked(tree *tr, int nodeA, int nodeB); 


#endif
