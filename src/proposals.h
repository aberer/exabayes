/**
   @file proposals.h

   @brief All proposals that make the chains of ExaBayes move in its
   space.
    
 */ 

#ifndef _PROPOSALS_H
#define _PROPOSALS_H

#include "nclConfigReader.h"


#ifdef __cplusplus
extern "C"{
#endif


void step(state *curstate); 
void normalizeProposalWeights(state *curstate);
void setupProposals(state *chain, initParamStruct *initParams); 

#ifdef __cplusplus
}
#endif


#endif
