/**
   @file  proposalStructs.h
   @brief Provides the main structs for ExaBayes. 

   This is still kind of ugly, will be redesigned soon. 
*/ 


#ifndef _BAYES_H
#define _BAYES_H


#include "proposalType.h"
#include "config.h"
#include "rng.h"
#include "tune.h"
#include "eval.h"
#include "stack.h"

/* TODO not enabled yet, was not such a good idea */
typedef struct
{
  proposal_type type;
  int wasAccepted;

  
  int model;

  union
  {
    double alpha;    
  } change;

} chainHistory ;



#define NUM_PROP_CATS 5 
typedef enum _cats {
  TOPOLOGY = 1  , 
  BRANCH_LENGTHS = 2, 
  FREQUENCIES = 3,
  SUBSTITUTION_RATES = 4  ,
  RATE_HETEROGENEITY = 5
} category_t; 



typedef struct
{
  int modelNum; 
  double alpha ;
 
  /* TODO only works with DNA */
  int numRates; 
  double substRates[6]; 

  int numFreqs; 
  double frequencies[4];

} perPartitionInfo; 		/* relevant info from pInfo  */



typedef struct 
{
  /* topology */
  topol *topo; 
  
  /* branch lengths  */
  double *branchLengths; 

  /* substitution parameter */
  perPartitionInfo *infoPerPart; 

} paramDump; 


typedef struct _pfun  proposalFunction; 




/* TODO will be replaced completely with paths  */
typedef struct 
{
  int whichBranch; /// helps us to remember which branch we manipulated.       
  double bls[NUM_BRANCHES]; /// branch lengths for saving 

  branch pruned; 		
  branch  insertBranch; /// the ids of both nodes  that originally constituted a branch before insertion 

  branch pruningBranch ; /// the branch prio to pruning => convertion: thisNode id is the direction we wandered to!
  double neighborBls[NUM_BRANCHES]; 
  double nextNeighborBls[NUM_BRANCHES]; 
} topoRecord; /// information that helps us to restore the previous state when doing topological moves or manipulating the branch lengths





typedef struct _state
{  
  tree * tr;
#if HAVE_PLL == 1   
  partitionList *partitions; 
#endif  
  proposalFunction *prevProposal;  /// only used for runtime efficiency. Is NULL, if we just saved/applied the state. 
  /* branch currentRoot; */ 

  boolean wasAccepted; 	/// for debug only  

  int id;   
  int couplingId;  /// indicates how hot the chain is (i = 0 => cold chain), may change!
  int currentGeneration;   
  
  double priorProb; /// the prior probability of the current state   => store in log form  

  
  lnlContainer lnl;  /// the current likelihood of the state  */

  double hastings;/// the proposal ratio 

  double penaltyFactor; //change to the probability of picking a proposal  

  proposalFunction **proposals; 
  int numProposals; 
  double *categoryWeights; 

  /* new stuff that we need when having multiple chains  */
  FILE *topologyFile; 
  FILE *outputParamFile; 

  /* RNG */
  randKey_t rKey;
  randCtr_t rCtr;

  /* saves the entire space in the parameter space  */
  paramDump dump;
} state;



struct _pfun 
{
  proposal_type ptype; 
  category_t category; 
  char *name; 

  double initWeight; 
  double currentWeight; 
  
  successCtr sCtr;  /// counts acceptance / rejection

  void (*apply_func)( state *chain, struct _pfun *pf ); /// modifies according to the proposal drawn
  void (*eval_lnl) (state *chain, struct _pfun *pf);  /// chooses the cheapest way to evaluate the likelihood  of a proposal 
  void (*reset_func)( state *chain, struct _pfun *pf );    /// only resets all cheap changes to the partition/tr strcuts => no evaluation (that's the goal at least)
  void (*autotune)(state *chain, struct _pfun *pf); 


  /**
     tunable parameters
  */
  union
  {
    double multiplier;
    double eSprStopProb; 
    double slidWinSize;  	
    double dirichletAlpha; 	/* TODO not used  */
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
  }remembrance; 

}; 


#endif
