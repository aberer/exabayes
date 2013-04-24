/**
   @file gibbsLikeTopo.h
   @brief Implements metropolized gibb-sampling-like topological moves. 
 */ 


#ifndef _GIBBS_LIKE_TOPO_H
#define _GIBBS_LIKE_TOPO_H

typedef struct  _insertWeight
{
  branch b; 
  double lnl; 
  double weightInFirst; 
  double weightInSecond; 
  boolean containedInFirst; 
  boolean containedInSecond; 
  double ratio; 		/* TODO ? */
  struct _insertWeight *next; 
} insertList ; 



void applyGuidedSPR(Chain *chain, proposalFunction *pf); 
void resetGuidedSPR(Chain *chain, proposalFunction *pf); 
void evalGuidedSPR(Chain *chain, proposalFunction *pf); 

#endif
