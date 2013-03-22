#ifndef _CHAIN_H
#define  _CHAIN_H

#include "nclConfigReader.h"




#define TOPO_RESTORE 0 
#define TOPO_SAVE 1 

void saveTreeStateToChain(state *chain); 
void applyChainStateToTree(state *chain); 
void traverseInitFixedBL(nodeptr p, int *count, tree *tr,  double z ); 
void initializeIndependentChains(tree *tr, analdef *adef, state **resultIndiChains); 
double getChainHeat(state *chain ); 

#endif
