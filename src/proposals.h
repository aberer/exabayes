#ifndef _PROPOSALS_H
#define _PROPOSALS_H

#include "nclConfigReader.h"

void step(state *curstate); 
void normalizeProposalWeights(state *curstate);
void resetSuccessCounters(state *chain); 
void setupProposals(state *chain, initParamStruct *initParams); 

#endif
