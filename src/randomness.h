/**
   @file randomness.h
   
   @brief Wraps random number generation functions. 

   @notice Do not generate random numbers otherwise! 

   notice new convention: 
   
   * chain-specific random numbers (e.g., proposals, acceptance) are
   created from randCtr_t and randKey_t and => this way we easily
   can reset the rng

   * random numbers that are needed in a global context (e.g.,
   switching chains) are created from a global rng. (see global
   variables)
     
   * for the ctr for the chain-specific stuf: first int is the
   generation, second a ctr (since we may need more numbers for each
   step, starting from 0)
*/ 


#ifndef _RANDOMNESS_H
#define _RANDOMNESS_H

#include "rng.h"
#include "branch.h"
#include "stack.h"


/* #ifdef __cplusplus */
/* extern "C"{ */
/* #endif */

/* todo replace that once */
typedef struct _rngState
{
  randKey_t key; 
  randCtr_t ctr; 
} rngState; 


void initLocalRng(state *theChain); 

int drawGlobalRandIntBound(int upperBound); 
randCtr_t drawGlobalRandInt();
double drawGlobalDouble01();
int drawRandInt(state *chain, int upperBound); 
double drawRandDouble01(state *chain);
double drawRandExp(state *chain, double lambda);
double drawRandBiUnif(state *chain, double x);
int drawSampleProportionally(state *chain,  double *weights, int numWeight ); 
void drawPermutation(state *chain, int* perm, int n); 
branch drawBranchUniform(state *chain); 
void generateRandomPath( state *chain ,stack *s, double stopProp); 
branch drawSubtreeUniform(state *chain); 
double drawMultiplier(state *chain, double multiplier); 
double drawFromSlidingWindow(state *chain, double param, double window); 
branch drawInnerBranchUniform(state *chain); 



/* #ifdef __cplusplus */
/* } */
/* #endif */

#endif
