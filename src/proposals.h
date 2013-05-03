/**
   @file proposals.h

   @brief All proposals that make the chains of ExaBayes move in its
   space.
    
 */ 

#ifndef _PROPOSALS_H
#define _PROPOSALS_H

#include "nclConfigReader.h"


void initProposalFunction( proposal_type type, initParamStruct *initParams, proposalFunction **result); 
double tuneParameter(int batch, double accRatio, double parameter, boolean inverse ); 
#endif
