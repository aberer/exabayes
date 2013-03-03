#ifndef _CHAIN_H
#define  _CHAIN_H

#include "nclConfigReader.h"




#define TOPO_RESTORE 0 
#define TOPO_SAVE 1 

void saveTreeStateToChain(state *chain, tree *tr); 
void applyChainStateToTree(state *chain, tree *tr); 
void traverseInitFixedBL(nodeptr p, int *count, tree *tr,  double z ); 
void initializeIndependentChains(tree *tr, state **resultIndiChains, initParamStruct **initParamsPtr); 
void printInfo(state *chain, const char *format, ...);
#endif
