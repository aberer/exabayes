/**
   @file  proposalStructs.h
   @brief Provides the main structs for ExaBayes. 

   This is still kind of ugly, will be redesigned soon. 
*/ 


#ifndef _PROPOSAL_STRUCTS_H
#define _PROPOSAL_STRUCTS_H


#include "proposalType.h"
#include "config.h"
#include "rng.h"




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





typedef struct
{
  int acc; 
  int rej;     
} accRejCtr; 


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
  
  /* to avoid that anything goes wrong */
  double likelihood; 
} paramDump; 


typedef struct _pfun  proposalFunction; 





typedef struct 
{
  int whichBranch; /// helps us to remember which branch we manipulated.       
  double bls[NUM_BRANCHES]; /// branch lengths for saving 

  int prunedSubTree ; 

  int insertBranch[2]; /// the ids of both nodes  that originally constituted a branch before insertion 

  int pruningBranches[2] ; /// the ids of the two neighbors prior to pruning 
  double neighborBls[NUM_BRANCHES]; 
  double nextNeighborBls[NUM_BRANCHES]; 
} topoRecord; /// information that helps us to restore the previous state when doing topological moves or manipulating the branch lengths





typedef struct _state
{
  tree * tr;
#if HAVE_PLL == 1   
  partitionList *partitions; 
#endif  
  
  double curprior;
  double newprior;
  double hastings;

  double penaltyFactor; //change to the probability of picking a proposal

  double likelihood; 


  proposalFunction **proposals; 
  int numProposals; 
  double *categoryWeights; 
  
  int currentGeneration; 

  /* indicates how hot the chain is (i = 0 => cold chain), may change! */
  int couplingId; 

  /* new stuff that we need when having multiple chains  */
  FILE *topologyFile; 
  FILE *outputParamFile; 
  /* unique id, also determines */
  int id; 

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
  
  accRejCtr successCtr; 

  void (*apply_func)( state *curstate, struct _pfun *pf );
  void (*reset_func)( state *curstate, struct _pfun *pf );   

  /* TODO not used yet */
  void (*autotune)(struct _pfun *pf); 

  /* TODO use something like that */
  /* double (*get_prior_ratio) (state *curstate);  */


  /**
     tunable parameters
  */
  union
  {
    double eSprStopProb; 
    double slidWinSize;  	
    double dirichletAlpha; 	/* TODO not used  */
    double stdDev ; 		/* TODO not used  */
  } parameters ; 
  /* more parameters could be added in another union. This is a bit of
     a hack for simplicity, but saves us all that memory allocation
     stuff */



  /**
     Variables that help us remember what we changed. 
   */
  union 
  {
    topoRecord *topoRec; 
    perPartitionInfo *partInfo; 
  }remembrance; 

}; 


#endif
