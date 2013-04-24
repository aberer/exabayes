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


/* todo replace that once */
typedef struct _rngState
{
  randKey_t key; 
  randCtr_t ctr; 
} rngState; 


void initLocalRng(Chain *theChain); 

int drawGlobalRandIntBound(int upperBound); 
randCtr_t drawGlobalRandInt();
double drawGlobalDouble01();
int drawRandInt(Chain *chain, int upperBound); 
double drawRandDouble01(Chain *chain);
double drawRandExp(Chain *chain, double lambda);
double drawRandBiUnif(Chain *chain, double x);
double drawRandGamma(Chain *chain, double alpha, double beta);
void drawRandDirichlet(Chain *chain, double* results, double* alphas, int length);
void drawDirichletExpected(Chain *chain, double* results, double* mean, double beta, int length);
int drawSampleProportionally(Chain *chain,  double *weights, int numWeight ); 
void drawPermutation(Chain *chain, int* perm, int n); 
branch drawBranchUniform(Chain *chain); 
void generateRandomPath( Chain *chain ,stack *s, double stopProp); 
branch drawSubtreeUniform(Chain *chain); 
double drawMultiplier(Chain *chain, double multiplier); 
double drawFromSlidingWindow(Chain *chain, double param, double window); 
branch drawInnerBranchUniform(Chain *chain); 


double densityDirichlet(double *values, double *alphas, int length);


#endif
