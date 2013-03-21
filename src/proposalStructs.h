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

/* #include "history.h" */

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







/* TODO sliding windows are not an attribute of remembrance */


/* spr proposal specific */
typedef struct {
  double qz[NUM_BRANCHES]; /* BL values prior to re-insertion */
  int single_bl_branch;
} branch_length_remembrance; 



/* model proposal specific  */
typedef struct
{
  analdef * adef;
  int model;
  int nstates;
  int numSubsRates;
  double *curSubsRates;//used for resetting
} model_remembrance; 

/* frequency proposal specific  */
typedef struct
{
  analdef * adef;
  int model;
  int numFrequRates;
  double *curFrequRates;//used for resetting
} frequency_remembrance; 


/* spr proposal specific  */
typedef struct 
{
  nodeptr p; /* node pruned */
  nodeptr nb;   /* p->next->back describes an edge that dissapears when p is pruned */
  double nbz[NUM_BRANCHES];
  nodeptr nnb; /* p->next->next->back describes an edge that dissapears when p is pruned */
  double nnbz[NUM_BRANCHES];
  nodeptr r; /* edge neighbour of re-insertion node, q->back */
  nodeptr q; /* re-insertion node */
} spr_move_remembrance; 



typedef struct
{
  double alpha ; 
  /* TODO only works with DNA */
  double substRates[6]; 
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


typedef struct _state
{
  /* TODO remove  */
  /* these 3 are independent of the state, can be taken out unless we want to pass a single pointer as an argument*/  
  nodeptr * list; /* list of possible re-insertion nodes */ 

  tree * tr;
#if HAVE_PLL == 1   
  partitionList *partitions; 
#endif  

  
  double curprior;
  double newprior;
  double hastings;

  double penaltyFactor; //change to the probability of picking a proposal

  double likelihood; 


  /* TODO that does not work currently */
  /* chainHistory history;  */


  /* TODO make all of this more generic in form of a history     */
  spr_move_remembrance sprMoveRemem; 
  branch_length_remembrance brLenRemem; 
  /* gamma_remembrance gammaRemem;  */
  model_remembrance modelRemem; 
  frequency_remembrance frequRemem; 


  proposalFunction **proposals; 
  int numProposals; 
  double *categoryWeights; 

  /* TODO should this be an attribute of the chain? we could have a run-struct instead...*/
  int numGen; 
  int currentGeneration; 
  /* int samplingFrequency;  */

  /* indicates how hot the chain is (i = 0 => cold chain), may change! */
  int couplingId; 

  /* new stuff that we need when having multiple chains  */
  FILE *topologyFile; 
  FILE *outputParamFile; 
  /* unique id, also determines */
  int id; 

  /* TODO 
     these things are needed for convergence diagnostics, but we
     should not store it in here

     this pointer is shared among all chains 
  */
  /* hashtable *bvHash;  */

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


  /* tunable parameters   */
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




  /* TODO wrap this up properly  */
  /*  not ideal yet */
  union 
  {
    double alpha; 
  }remembrance; 
  int model; 




}; 


#endif
