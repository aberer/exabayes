/**
   @file  proposalStructs.h
   @brief Provides the main structs for ExaBayes. 

   This is still kind of ugly, will be redesigned soon. 
*/ 


#ifndef _BAYES_H
#define _BAYES_H


using namespace std; 
#include <iostream>

#include "proposalType.h"
#include "config.h"
#include "rng.h"
#include "eval.h"
#include "stack.h"
#include "SuccessCtr.hpp"
#include "categoryType.h"
#include "Chain.hpp"
#include "nclConfigReader.h"


class TreeAln; 
class LnlRestorer; 

typedef struct _pfun  proposalFunction; 


/* TODO will be replaced completely with paths  */
typedef struct 
{
  /* int whichBranch; /// helps us to remember which branch we manipulated.        */
  double bls[NUM_BRANCHES]; /// branch lengths for saving 

  branch pruned; 		
  branch  insertBranch; /// the ids of both nodes  that originally constituted a branch before insertion 

  branch pruningBranch ; /// the branch prio to pruning => convertion: thisNode id is the direction we wandered to!
  double neighborBls[NUM_BRANCHES]; 
  double nextNeighborBls[NUM_BRANCHES]; 
} topoRecord; /// information that helps us to restore the previous Chain when doing topological moves or manipulating the branch lengths


struct _pfun 
{
  proposal_type ptype; 
  category_t category; 
  char *name; 

  double initWeight; 
  double currentWeight; 

  void (*apply_func)( Chain *chain, struct _pfun *pf ); /// modifies according to the proposal drawn
  void (*eval_lnl) (Chain *chain, struct _pfun *pf);  /// chooses the cheapest way to evaluate the likelihood  of a proposal 
  void (*reset_func)( Chain *chain, struct _pfun *pf );    /// only resets all cheap changes to the partition/tr strcuts => no evaluation (that's the goal at least)
  void (*autotune)(struct _pfun *pf, SuccessCtr *ctr); 



  double relativeWeight; 

  /**
     tunable parameters
  */
  union
  {
    double multiplier;
    double eSprStopProb; 
    double slidWinSize;  	
    int radius;     		/* for guided spr moves */
    double dirichletAlpha;
    double stdDev ; 		/* TODO not used  */
  } parameters ; 
  /* more parameters could be added in another union. This is a bit of
     a hack for simplicity, but saves us all that memory allocation
     stuff */


  union
  {
    double multiplier ; 	/* the lambda of eSPR for multiplying moves  */
  } param2;  

  /**
     Variables that help us remember what we changed. 
   */
  union 
  {
    stack *modifiedPath; 
    topoRecord *topoRec;
    perPartitionInfo *partInfo; 
    double multiplier; 
  }remembrance;

  /* TODO dirty: is also a remembrance variable  */
  double ratio; 
}; 

void exa_main(analdef *adef, int seed, initParamStruct *initParam); 
void setupGlobals(initParamStruct *initParams); 

#endif
