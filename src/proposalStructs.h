#ifndef _PROPOSAL_STRUCTS_H
#define _PROPOSAL_STRUCTS_H



/* okay, so defining enums this way is rather save  */
#define NUM_PROPOSALS (6)
typedef enum
{
E_SPR = 0,
UPDATE_MODEL = 1 ,
UPDATE_GAMMA = 2 ,
UPDATE_GAMMA_EXP=3,
UPDATE_SINGLE_BL = 4,
UPDATE_SINGLE_BL_EXP = 5
} proposal_type;




/* TODO sliding windows are not an attribute of remembrance */


/* spr proposal specific */
typedef struct {
  double qz[NUM_BRANCHES]; /* BL values prior to re-insertion */
  double bl_sliding_window_w;
  double bl_prior;
  double bl_prior_exp_lambda;
  int single_bl_branch;  
} branch_length_remembrance; 


/* gamma proposal specific  */
typedef struct
{
  double curAlpha;
  double gm_sliding_window_w;  
} gamma_remembrance; 


/* model proposal specific  */
typedef struct
{
  analdef * adef;
  int model;
  int nstates;
  int numSubsRates;
  double rt_sliding_window_w;
  double *curSubsRates;//used for resetting
} model_remembrance; 


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
  /* these 3 are independent of the state, can be taken out unless we want to pass a single pointer as an argument*/  
  nodeptr * list; /* list of possible re-insertion nodes */ 
  int maxradius;  /* maximum radius of re-insertion from the pruning point */
  tree * tr;
  
  double curprior;
  double newprior;
  double hastings;
  
  double penaltyFactor; //change to the probability of picking a proposal

  double eSprStopProb; 

  spr_move_remembrance sprMoveRemem; 
  branch_length_remembrance brLenRemem; 
  gamma_remembrance gammaRemem; 
  model_remembrance modelRemem; 

  double proposalWeights[NUM_PROPOSALS]; 
  double proposalLogisticT[NUM_PROPOSALS];
  /* prob_bucket_t proposals[NUM_PROPOSALS];  */
  int acceptedProposals[NUM_PROPOSALS]; 
  int rejectedProposals[NUM_PROPOSALS];
  int totalAccepted;
  int totalRejected;

  /* TODO should this be an attribute of the chain? we could have a run-struct instead...*/
  int numGen; 
  int currentGeneration; 
  int samplingFrequency; 
} state;




typedef struct {
  int  ptype;  
  void (*apply_func)( state * curstate );
  void (*reset_func)( state * curstate );    
  double (*prior_function) (void *something); 

} proposal_functions;


#endif

