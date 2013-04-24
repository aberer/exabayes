
/**
   @file chain.h
   @brief Functions that are applied to individual chains. 
*/ 


#ifndef _CHAIN_OLD_H
#define _CHAIN_OLD_H

#include "nclConfigReader.h"

#define TOPO_RESTORE 0 
#define TOPO_SAVE 1 

class Chain; 
class TreeAln; 
class CoupledChains; 

void saveTreeStateToChain(Chain *chain); 
void applyChainStateToTree(Chain *chain); 
void traverseInitFixedBL(nodeptr p, int *count, TreeAln *traln,  double z ); 
void initializeIndependentChains( analdef *adef, int seed, vector<CoupledChains*> &resultIndiChains, initParamStruct *initParams ); 

void drawProposalFunction(Chain *chain, proposalFunction **result ); 

void initParamDump(TreeAln *traln, paramDump *dmp); 
void copyParamDump(TreeAln *traln,paramDump *dest, const paramDump *src); 

#endif
