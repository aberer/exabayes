#include "axml.h"
#include "GlobalVariables.hpp"
#include "output.h"
#include "eval.h"
#include "adapters.h"
#include "branch.h"
#include "Path.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "densities.h"

void expensiveVerify(tree *tr); 

void edit_subs_rates(Chain *chain, int model, int subRatePos, double subRateValue);


#if 0 

/* quarantining this, until change is over */

//reads proposalWeights p and sets t for logistic function such that f(t)=p
void findLogisticT(Chain *chain){
  static double leftShift=2.0;
   for(int i = 0; i < NUM_PROPOSALS;++i){
     if(chain->proposalWeights[i]>0){
    chain->proposalLogisticT[i]=-log((1.0-chain->proposalWeights[i])/chain->proposalWeights[i])+leftShift;  
    //printf("p(i): %f e^..: %f\n",chain->proposalWeights[i], (1.0/(1.0+exp(leftShift-chain->proposalLogisticT[i])) ));
     }else{
       chain->proposalLogisticT[i]=-1.0;
     }
   }
}

//get values for logistic function
void findLogisticP(Chain *chain){
  static double leftShift=2.0;
  for(int i = 0; i < NUM_PROPOSALS;++i){
    if(chain->proposalLogisticT[i]>=0){
    chain->proposalWeights[i]=(1.0/(1.0+exp(leftShift-chain->proposalLogisticT[i])) );  
  //  printf("p(i): %f e^..: %f t: %f\n",chain->proposalWeights[i], (1.0/(1.0+exp(leftShift-chain->proposalLogisticT[i])) ), chain->proposalLogisticT[i]);
   }
  }
}
#endif


static void recordSubsRates(Chain *chain, int model, int numSubsRates, double *prevSubsRates)
{
  pInfo *partition = chain->traln->getPartition(model); 
  assert(partition->dataType == DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = chain->traln->getSubstRate(model,i); 
}


static void record_branch_info(TreeAln* traln, nodeptr p, double *bl)
{
  int i;
  int numBranches = traln->getNumBranches(); 
  for(i = 0; i < numBranches; i++)
    bl[i] = traln->getBranchLength( p,0);
}


static void onePartitionEval(Chain *chain, proposalFunction *thisProposal)
{
  int model = thisProposal->remembrance.partInfo->modelNum;
  branch root = findRoot(chain->traln->getTr());  
  nodeptr p = findNodeFromBranch(chain->traln->getTr(), root); 
  evaluateOnePartition(chain, p, TRUE, model); 
}


static void evalBranch(Chain *chain, proposalFunction *thisProposal)
{
  assert(thisProposal->remembrance.modifiedPath->size() == 1); 
  branch b = thisProposal->remembrance.modifiedPath->at(0);
  nodeptr p = findNodeFromBranch(chain->traln->getTr(), b); 
  evaluateGenericWrapper(chain,p, FALSE ); 
}


static void dummy_eval(Chain *chain,proposalFunction *thisProposal)
{
  evaluateGenericWrapper(chain, chain->traln->getTr()->start, TRUE );
}


/* is this function really uniform? a tip node has a lower
   propbability to be drawn (only one path leads to it) */ 
static nodeptr select_random_subtree(Chain *chain, tree *tr)
{
  nodeptr 
    p;

  do
    {
      int 
        exitDirection = chain->getChainRand()->drawRandInt(3); 
     
      int r = chain->getChainRand()->drawRandInt(tr->mxtips - 2) ; 

      p = tr->nodep[ r + 1 + tr->mxtips];
      
      switch(exitDirection)
        {
        case 0:
          break;
        case 1:
          p = p->next;
          break;
        case 2:
          p = p->next->next;
          break;
        default:
          assert(0);
        }
    }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));

  assert(!isTip(p->number, tr->mxtips));

  return p;
}


//--------Alpha-Proposal-for-GAMMA-----------------------------------------------------------------

double get_alpha_prior(Chain *chain )
{  
 return 1;//TODO obviously needs acctual prior 
}


static void simple_gamma_proposal_apply(Chain * chain, proposalFunction *pf)
{
  //TODO: add safety to max and min values
  double newalpha, curv, r,mx,mn;

  int model = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions());  
  perPartitionInfo* remem = pf->remembrance.partInfo; 
  remem->modelNum = model; 

  curv = chain->traln->getAlpha(model); 
  remem->alpha = curv; 

  switch(pf->ptype)
    {
    case UPDATE_GAMMA:
      {
	double slidWin = pf->parameters.slidWinSize;
	    
	/* case STANDARD://simple sliding window */
	r = chain->getChainRand()->drawRandDouble01();
	mn = curv-(slidWin/2);
	mx = curv+(slidWin/ 2);
	newalpha = fabs(mn + r * (mx-mn));
      }
      break;
    case UPDATE_GAMMA_EXP: 
      newalpha  = chain->getChainRand()->drawRandExp(1/curv);
      break;
    default:
      assert(0);
    }
  
  chain->hastings = 1; //since it is symmetrical, hastings=1
  
  chain->traln->setAlphaSave(newalpha, model); 
  chain->traln->discretizeGamma(model); 
}


static void simple_gamma_proposal_reset(Chain *chain, proposalFunction *pf)
{
  perPartitionInfo *info  = pf->remembrance.partInfo;     
  chain->traln->setAlphaSave(info->alpha, info->modelNum); 
  chain->traln->discretizeGamma(info->modelNum); 
}

//------------------------------------------------------------------------------


#if 0 
/* TODO we do not even use this function, do we? NOTE Now we do ;) */
void penalize(Chain *chain, int which_proposal, int acceptance)
{
  double max=4.0;
  double min=0.0;
  if(acceptance){
    chain->proposalLogisticT[which_proposal]+=chain->penaltyFactor;
    if(chain->proposalLogisticT[which_proposal]>max)
      chain->proposalLogisticT[which_proposal]=max;
    //chain->proposalWeights[which_proposal] /= chain->penaltyFactor; 
  }
  else {
  /*  if(chain->totalRejected>0){
    chain->proposalLogisticT[which_proposal]-=((chain->totalAccepted/chain->totalRejected)*chain->penaltyFactor);*///TODO check whether "VCG" is better than considering all
    if(chain->totalRejected - chain->rejectedProposals[which_proposal]>0){
    chain->proposalLogisticT[which_proposal]-=(((chain->totalAccepted- chain->acceptedProposals[which_proposal])/(chain->totalRejected - chain->rejectedProposals[which_proposal]))*chain->penaltyFactor);      
    }else{
      chain->proposalLogisticT[which_proposal]-=chain->penaltyFactor;
    }
    if(chain->proposalLogisticT[which_proposal]<min)
      chain->proposalLogisticT[which_proposal]=min;
    //chain->proposalWeights[which_proposal] *= chain->penaltyFactor; 
  }
  findLogisticP(chain);
  //normalizeProposalWeights(chain); 
}

#endif


static void simple_model_proposal_apply(Chain *chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 

  int model = chain->getChainRand()->drawRandInt( chain->traln->getNumberOfPartitions());
  perPartitionInfo *info = pf->remembrance.partInfo;

  pInfo *partition = chain->traln->getPartition( model); 
  int numRates = (partition->states * partition->states - partition->states) / 2; 

  recordSubsRates(chain, model , numRates, info->substRates);

  info->modelNum = model;
  info->numRates = numRates; 

  //choose a random set of model params,
  //probably with dirichlet proposal
  //with uniform probabilities, no need to have other
  int state, changeState = 0;
  double new_value = 0,curv;
  double r;
  
  double mx,mn; //for sliding window
  
  int* list  = (int*)exa_calloc(numRates, sizeof(int)); 
  // int list[numRates];//for biunif_distr and biunif_perm_distr
  
  /* int numberOfEdits;//for biunif_perm_distr */
  
  proposal_type  pType = pf->ptype; 


  /* TODO reactivate when used and replace old types   */
  /* switch(pType) */
  /*   { */
  /*   case STANDARD: */
  /*   case BIUNIF_DISTR: */
  /*     /\* numberOfEdits=chain->modelRemem.numSubsRates; *\/ */
  /*     break; */
  /*   case BIUNIF_PERM_DISTR: */
  /*     /\* numberOfEdits=drawRandInt(chain, chain->modelRemem.numSubsRates); *\/ */
  /*     break; */
  /*   default : assert(0);  */
  /*   }  */

  for(state = 0;state<numRates ; state ++)
    {      
      switch(pType)
        {
	case UPDATE_MODEL : //using the branch length sliding window for a test    
	  {
	    double range = pf->parameters.slidWinSize / 2; 
	    
	    changeState=state;
	    curv = traln->getSubstRate(model, state);
	    r =  chain->getChainRand()->drawRandDouble01();
	    mn = curv-range ;
	    mx = curv+range ;
	
	    new_value = fabs(mn + r * (mx-mn));
	  
	    /* Ensure always you stay within this range */
	    if(new_value > RATE_MAX) new_value = RATE_MAX;
	    if(new_value < RATE_MIN) new_value = RATE_MIN;
	  }
	  break;
        
	case UPDATE_MODEL_BIUNIF:
	  changeState=chain->getChainRand()->drawRandInt( numRates);
	  if(list[changeState]!=1)
	    {
	      list[changeState]=1;;      
	      curv = traln->getSubstRate(model, changeState); 
	      r =  chain->getChainRand()->drawRandBiUnif( curv);
	      new_value = r;
	      chain->hastings*=curv/new_value;
 
	      while(new_value> RATE_MAX|| new_value< RATE_MIN)
		{
		  if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
		  if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
		}
	    }
	  break;
	  
	  

	  /* TODO activate and replace */
	/* case BIUNIF_PERM_DISTR://TODO NOT used. Lower values than before (without subType) */
	/*   drawPermutation(chain, list, chain->modelRemem.numSubsRates); */
	/*   curv = partition->substRates[list[state]]; */
	/*   r =  drawRandBiUnif(chain, curv); */
	/*   new_value = r; */
	/*   while(new_value> RATE_MAX|| new_value< RATE_MIN){ */
	/*     if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value; */
	/*     if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value; */
	/*   } */
      
	/*   chain->hastings*=curv/new_value; */
	/*   break; */
	 
	/* case SINGLE_BIUNIF://TODO not used. Figure out error */
	/*   drawPermutation(chain,list, chain->modelRemem.numSubsRates); */
	/*   curv = partition->substRates[list[state]]; */
	/*   r =  drawRandBiUnif(chain, curv); */
	/*   new_value = r; */
	/*   while(new_value> RATE_MAX|| new_value< RATE_MIN){ */
	/*     if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value; */
	/*     if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value; */
	/*   } */
      
	/*   chain->hastings*=curv/new_value; */
	/*   break; */
	 
        default:
	  {
	    assert(0);
	  }

	}//end of switch
      edit_subs_rates(chain, info->modelNum, changeState, new_value);
    }
  //recalculate eigens

  chain->traln->initRevMat(model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */


  exa_free(list); 
    
  /* TODO: need to broadcast rates here for parallel version ! */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/ */

  //TODO: without this, the run will fail after a successful model, but failing SPR
  //TODOFER: what did we have in mind regarding the comment above?
  
  /* evaluateGeneric(chain->tr, chain->tr->start, FALSE); */
  //for prior, just use dirichlet
  // independent gamma distribution for each parameter
  //the pdf for this is
  // for gamma the prior is gamma

  //for statefreqs should all be uniform

  //only calculate the new ones
}



static void model_dirichlet_proposal_apply(Chain *chain, proposalFunction *pf)//llpqr
{

  int model = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions());
  perPartitionInfo *info = pf->remembrance.partInfo;

  pInfo *partition = chain->traln->getPartition( model); 
  int numRates = (partition->states * partition->states - partition->states) / 2; 

  recordSubsRates(chain, model , numRates, info->substRates);

  info->modelNum = model;
  info->numRates = numRates; 

  double r[numRates];

  int* list  = (int*)exa_calloc(numRates, sizeof(int)); 

     
  chain->getChainRand()->drawDirichletExpected(r, partition->substRates, pf->parameters.dirichletAlpha, numRates);
  chain->hastings= densityDirichlet(partition->substRates, r, numRates) 
    / densityDirichlet(r, partition->substRates, numRates);    
  
  /* Ensure always you stay within this range */
  for(int i=0; i<numRates; i++)
    {
      if(r[i] > RATE_MAX) r[i] = RATE_MAX;
      if(r[i] < RATE_MIN) r[i] = RATE_MIN;
      edit_subs_rates(chain, info->modelNum, i, r[i]);
    }
     
  //recalculate eigens

  chain->traln->initRevMat(model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */


  exa_free(list); 

}


// //draws a random subset (drawing with replacement) of the states and changes the according to biunif distribution.
// static void biunif_model_proposal_apply(state *chain, int pSubType)
// {
//   //record the old one 
//    chain->modelRemem.model=drawRandInt(chain->tr->NumberOfModels);
//   recordSubsRates(chain->tr, chain->modelRemem.model, chain->modelRemem.numSubsRates, chain->modelRemem.curSubsRates);
//   int state, randState;
//   double new_value,curv;
//   double r;
//   int list[chain->modelRemem.numSubsRates];
//   
//   
//    chain->hastings=1.0;
//   for(state = 0;state<chain->modelRemem.numSubsRates ; state ++)
//     {
//       randState=drawRandInt(chain->modelRemem.numSubsRates);
//       if(list[randState]!=1)
//       {
//       list[randState]=1;;
//       
//       curv = chain->tr->partitionData[chain->modelRemem.model].substRates[randState];
//       r =  drawRandBiUnif(curv);
// 
//       new_value = r;
// 
//       while(new_value> RATE_MAX|| new_value< RATE_MIN){
//       if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
//       if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
//       }
//       
//        chain->hastings*=curv/new_value;
//       edit_subs_rates(chain->tr,chain->modelRemem.model, randState, new_value);
//       }
//     }
//       
// 
// 
//   initReversibleGTR(chain->tr, chain->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
// 
//   evaluateGeneric(chain->tr, chain->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
// 
// }



static void perm_biunif_model_proposal_apply(Chain *chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 

  int model = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions());; 

  pInfo *partition = chain->traln->getPartition( model) ; 

  perPartitionInfo *info = pf->remembrance.partInfo; 

  int numRates  = (partition->states  * partition->states -partition->states ) / 2 ;  
  info->numRates = numRates; 

  recordSubsRates(chain, model, numRates, info->substRates);
  int state, randNumber;
  double new_value,curv;
  double r;
  
  randNumber=chain->getChainRand()->drawRandInt(numRates);
  int *perm = (int*)exa_calloc(numRates, sizeof(int)); 
  // int perm[numRates];
  chain->getChainRand()->drawPermutation(perm, numRates);

  for(state = 0;state < randNumber ; state ++)
    {           
      curv = traln->getSubstRate(model, perm[state]);;
      r =  chain->getChainRand()->drawRandBiUnif(curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
	if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
	if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
      chain->hastings*=curv/new_value;
      edit_subs_rates(chain,model, perm[state], new_value);      
    }
      
  chain->traln->initRevMat(model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/   */
  exa_free(perm); 
}



static void single_biunif_model_proposal_apply(Chain *chain,proposalFunction *pf)//NOTE whenever a model parameter changes, all branch lengths have to be re-normalized with 1/fracchange. Additionally we always must do a full tree traversal to get the likelihood. So updating a single parameter is rather expensive, .
{
  TreeAln *traln = chain->traln; 

  //record the old one //TODO sufficient to store single value.  
  int model = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions()); 

  perPartitionInfo *info =  pf->remembrance.partInfo; 
  
  pInfo *partition = chain->traln->getPartition(model) ; 

  info->modelNum = model; 
  info->numRates = partition->states ;   
  int numRates = info->numRates; 
  recordSubsRates(chain, model, info->numRates, info->substRates);
  //choose a random set parameter,
  //with uniform probabilities

  int randState= chain->getChainRand()->drawRandInt(numRates);

  double new_value,curv;
  double r;
  
  //int state=drawRandInt(chain->modelRemem.numSubsRates);
  

  curv = traln->getSubstRate(model, randState);
  r =  chain->getChainRand()->drawRandBiUnif(curv);

  new_value = r;
      
  while(new_value> RATE_MAX|| new_value< RATE_MIN){
    if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
    if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
  }

  edit_subs_rates(chain,model, randState, new_value);

  chain->hastings=curv/new_value;

  
  chain->traln->initRevMat(model);
  // exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
  
  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/ */
}

static void all_biunif_model_proposal_apply(Chain *chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 
  
  int model = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions());

  perPartitionInfo *info = (perPartitionInfo* ) pf->remembrance.partInfo; 
  
  info->modelNum = model; 
  
  //record the old one 
  pInfo *partition = chain->traln->getPartition(model) ; 
  info->numRates = (partition->states * partition->states - partition->states   ) / 2 ; 
  int numRates = info->numRates; 

  recordSubsRates(chain, model, numRates, info->substRates);
  //choose a random set parameter,
  //with uniform probabilities
  int state;
  double new_value,curv;
  double r;

  for(state = 0;state<numRates ; state ++)
    {
      curv = traln->getSubstRate(model, state); 
      r =  chain->getChainRand()->drawRandBiUnif(curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       chain->hastings*=curv/new_value;
      edit_subs_rates(chain,model, state, new_value);
    }
  chain->traln->initRevMat(model);
  // exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/ */
}






static void restore_subs_rates(Chain *chain, int model, int numSubsRates, double *prevSubsRates)
{
  TreeAln *traln = chain->traln; 

  pInfo *partition = chain->traln->getPartition( model); 

  assert(partition->dataType == DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)	
    traln->setSubstSave( prevSubsRates[i],model,i ); 

  chain->traln->initRevMat(model);
  // exa_initReversibleGTR(chain, model);

  /* TODO need to broadcast rates here for parallel version */

  /* evaluateOnePartition(chain, tr->start, TRUE, model);  */
}


//--------Branch-Length_Proposals---------------------------------------------------

double get_branch_length_prior( Chain *chain)
{//TODO decide on sensible prior
  return 1;  
}

//setting this out to allow for other types of setting
static void set_branch_length_sliding_window(Chain *chain, nodeptr p, int numBranches,Chain * s, boolean record_tmp_bl, double windowRange, double *bls)
{
  TreeAln *traln = chain->traln; 
  int i;
  double newZValue;
  double r,mx,mn;
  /* double newValue = 0; */

  for(i = 0; i < numBranches; i++)
    {
      double real_z;
    
      if(record_tmp_bl)
	{
	  assert(p->z[i] == p->back->z[i]);	  	  
	  bls[i] = traln->getBranchLength( p,0); 

	  /* bls[i] = p->back->z_tmp[i] = p->z[i]; */
	  /* p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /\* keep current value *\/ */
	}
      r = chain->getChainRand()->drawRandDouble01();

      /* TODO@kassian just wondering: is it correct to use the global
	 frac-change here? what are the local frac-changes good for
	 then?  */
      real_z = -log( traln->getBranchLength( p,0)) * chain->traln->getTr()->fracchange;
      
      //     printf( "z: %f %f\n", p->z[i], real_z );
    
      mn = real_z- windowRange;
      mx = real_z+ windowRange;
      newZValue = exp(-(fabs(mn + r * (mx-mn)/chain->traln->getTr()->fracchange )));

      /* newValue = mn + r * (mx-mn); */
      /* s->newprior=get_branch_length_prior(&newValue); */
    
      /* Ensure always you stay within this range */
      if(newZValue > zmax) newZValue = zmax;
      if(newZValue < zmin) newZValue = zmin;
      assert(newZValue <= zmax && newZValue >= zmin);
      //     printf( "z: %f %f %f\n", p->z[i], newZValue, real_z );
      traln->setBranchLengthSave(newZValue,i,p); 
      // p->z[i] = p->back->z[i] = newZValue;
      //assuming this will be visiting each node, and multiple threads won't be accessing this
      //s->bl_prior += log(exp_pdf(s->bl_prior_exp_lambda,newZValue));
      //s->bl_prior += 1;
    }
}

static void set_branch_length_biunif(Chain *chain, nodeptr p, int numBranches,Chain * s, boolean record_tmp_bl, double *bls)
{
  TreeAln *traln = chain->traln; 
  int i;
  double newZValue;
  double r;
  double newValue = 0; 

  for(i = 0; i < numBranches; i++)
    {
      double real_z;
    
      if(record_tmp_bl)
	{
	  assert(p->z[i] == p->back->z[i]); 
	  bls[i] = traln->getBranchLength( p,0); 
	  /* p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /\* keep current value *\/ */
	}
      
      real_z = -log( traln->getBranchLength( p,0)) * chain->traln->getTr()->fracchange; //convert from exponential to real form
      r = chain->getChainRand()->drawRandBiUnif( real_z);//draw from [real_z/2,2*real_z]

    
      newZValue = exp(-(fabs( r/chain->traln->getTr()->fracchange )));//convert to exponential form
      
      newValue = r; 
      
      s->hastings=real_z/newValue;
      /* s->newprior=get_branch_length_prior(&newValue); */
    
      /* Ensure always you stay within this range */
      /*
      if(newZValue > zmax) newZValue = zmax;
      if(newZValue < zmin) newZValue = zmin;
      assert(newZValue <= zmax && newZValue >= zmin);

      
      p->z[i] = p->back->z[i] = newZValue;
      */
      
     /* Ensure always you stay within this range, reflect if necessary */
      while(newZValue > zmax || newZValue < zmin){
      if(newZValue > zmax) newZValue = 2*zmax-newZValue;
      if(newZValue < zmin) newZValue = 2*zmin-newZValue;
      }
      assert(newZValue <= zmax && newZValue >= zmin);
      traln->setBranchLengthSave(newZValue,i,p); 

    }
}

static void set_branch_length_exp(Chain *chain, nodeptr p, int numBranches,Chain * s, boolean record_tmp_bl, double *bls)
{
  TreeAln *traln = chain->traln;
  int i;
  double newZValue;
  double r;
  double real_z;
  for(i = 0; i < numBranches; i++)
    {
     // double real_z;
    
      if(record_tmp_bl)
	{
	  assert(p->z[i] == p->back->z[i]); 
	  bls[i] = traln->getBranchLength( p,0); 
	  /* p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /\* keep current value *\/ */
	}
   //   r = drawRandExp(lambda);
      real_z = -log( traln->getBranchLength( p,0)) * chain->traln->getTr()->fracchange;
       r = chain->getChainRand()->drawRandExp( 1.0/real_z);

	s->hastings=(1/r)*exp(-(1/r)*real_z)/((1/real_z)*exp(-(1/real_z)*r));

	/* s->newprior=get_branch_length_prior(&r); */
	
	newZValue = exp(-(fabs(r /chain->traln->getTr()->fracchange )));
   
    
      /* Ensure always you stay within this range, reflect if necessary *///TODO reflects transformed bl, needs to reflect actual bl
      while(newZValue > zmax || newZValue < zmin){
      if(newZValue > zmax) newZValue = 2*zmax-newZValue;
      if(newZValue < zmin) newZValue = 2*zmin-newZValue;
      }
      assert(newZValue <= zmax && newZValue >= zmin);
      traln->setBranchLengthSave(newZValue, i,p); 
				// p->z[i] = p->back->z[i] = newZValue;
     }
     
}


static node *select_branch_by_id_dfs_rec( node *p, int *cur_count, int target, Chain *s ) {
  if( (*cur_count) == target ) 
    {
      return p;
    }
  
  if( !isTip( p->number, s->traln->getTr()->mxtips )) {
    // node *q = p->back;
    node *ret = NULL;
    
    ++(*cur_count);
    ret = select_branch_by_id_dfs_rec( p->next->back, cur_count, target, s );
    if( ret != NULL ) {
      return ret;
    }
    
    ++(*cur_count);
    ret = select_branch_by_id_dfs_rec( p->next->next->back, cur_count, target, s );
    
    return ret;
  } else {
    return NULL;
  }
}



static node *select_branch_by_id_dfs( node *p, int target, Chain *s ) {
  const int num_branches = (2 * s->traln->getTr()->mxtips) - 3;
  
  //   if( target == num_branches ) {
  //     return NULL;
  //   } 
  
  assert( target < num_branches );
  
  int cur_count = 0;
  node *ret = NULL;
  
  
  ret = select_branch_by_id_dfs_rec( p, &cur_count, target, s );
  
  if( ret != NULL ) {
    return ret;
  } else {
    // not in subtree below p. search on in the subtrees below p->back
  
    node *q = p->back;
    
    if( isTip( q->number, s->traln->getTr()->mxtips ) ) {
      assert( 0 ); // error: target must be invalid 'branch id'
      return NULL; 
    }
    
    ++cur_count;
    ret = select_branch_by_id_dfs_rec( q->next->back, &cur_count, target, s );
  
    if( ret != NULL ) {
      return ret;
    }
    
    
    ++cur_count;
    ret = select_branch_by_id_dfs_rec( q->next->next->back, &cur_count, target, s );
  
    return ret;
  
  }  
  assert(0);    
}




static void random_branch_length_proposal_apply(Chain * chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 

  int multiBranches = chain->traln->getNumBranches(); 
  assert(multiBranches == 1 ); 
  int target_branch = chain->getChainRand()->drawRandInt((chain->traln->getTr()->mxtips * 2) - 3); 

  node *p = select_branch_by_id_dfs( chain->traln->getTr()->start, target_branch, chain );
  double z = traln->getBranchLength( p,0); 
  double oldZ = 0; 
   
  
  switch(pf->ptype)
    {
    case UPDATE_SINGLE_BL: 
      //for each branch get the current branch length
      //pull a uniform like
      //x = current, w =window
      //uniform(x-w/2,x+w/2)
      set_branch_length_sliding_window(chain,p, multiBranches, chain, TRUE, pf->parameters.slidWinSize, &oldZ); 
      break;
	
    case UPDATE_SINGLE_BL_EXP: 
      set_branch_length_exp(chain,p, multiBranches, chain, TRUE,&oldZ);
      break;
    case UPDATE_SINGLE_BL_BIUNIF: 

      set_branch_length_biunif(chain, p, multiBranches, chain, TRUE,&oldZ);
      break;
    default:
      assert(0);
    }

  /* clearStack(pf->remembrance.modifiedPath);  */
  /* assert(pf->remembrance.modifiedPath != NULL);  */
  branch b = constructBranch(p->number, p->back->number); 
  pf->remembrance.modifiedPath->append( b); 
  
  branch &something = pf->remembrance.modifiedPath->at(0); 
  something.length[0] = z; 
}


// static void biunif_branch_length_proposal_apply(state * chain, int pSubType)
// {
//    
//   //for one branch get the current branch length
//   //pull a uniform like
//   //x = current,
//   //uniform(x/2,x*2)
//   
//   
//   const int num_branches = (chain->tr->mxtips * 2) - 3;
//   int target_branch = drawRandInt(num_branches); 
//   node *p = select_branch_by_id_dfs( chain->tr->start, target_branch, chain );
//  
//   //set_branch_length_sliding_window(p, chain->tr->numBranches, chain, TRUE);
//   set_branch_length_biunif(p, chain->tr->numBranches, chain, TRUE);
// 
//   chain->brLenRemem.single_bl_branch = target_branch;
//   evaluateGeneric(chain->tr, chain->tr->start, TRUE); /* update the tr->likelihood *///TODO see below
//   //   return TRUE;
// }

// static void exp_branch_length_proposal_apply(state * chain, int pSubType)
// {
//   const int num_branches = (chain->tr->mxtips * 2) - 3;
//   int target_branch = drawRandInt(num_branches); 
//   node *p = select_branch_by_id_dfs( chain->tr->start, target_branch, chain );
//   
//   set_branch_length_exp(p, chain->tr->numBranches, chain, TRUE);
// 
//   chain->brLenRemem.single_bl_branch = target_branch;
//   evaluateGeneric(chain->tr, chain->tr->start, TRUE); /* update the tr->likelihood *///TODO see below
// }

static void random_branch_length_proposal_reset(Chain * chain, proposalFunction *pf)
{
  /* branch b = pf->remembrance.modifiedPath->content[0];  */
  branch b = pf->remembrance.modifiedPath->at(0); 
  nodeptr p = findNodeFromBranch(chain->traln->getTr(), b); 
  chain->traln->setBranchLengthSave(b.length[0],0,p); 
  pf->remembrance.modifiedPath->clear(); 

  // ok, maybe it would be smarter to store the node ptr for rollback rather than re-search it...
  /* p = select_branch_by_id_dfs( chain->tr->start, pf->remembrance.topoRec->whichBranch, chain ); */  
  /* reset_branch_length(p, getNumBranches(chain->tr), pf->remembrance.topoRec->bls); */
  //   printf( "reset bl: %p %f\n", p, p->z[0] );
  //update_all_branches(chain, TRUE);

  /* i am not so sure, if we really can do the FALSE here ; will raxml recognize that a branch length has changed? disabling it for now... 

     TODO I think, we should evaluate at the respctive node 
   */

 // evaluateGeneric(chain->tr, p, FALSE); //This yields a very slight likelihood difference.NOTE if we want exact likelihoods as before the proposal, we must evaluate from chain->tr->start, that is: evaluateGeneric(chain->tr, chain->tr->start, TRUE);
}

double get_frequency_prior(Chain * chain)
{
 return 1; 
}

static void restore_frequ_rates(Chain *chain, int model, int numFrequRates, double *prevFrequRates)
{
  /* NOTICE: this function should not be called repeatedly  */

  pInfo *partition = chain->traln->getPartition(model);

  assert(partition->dataType == DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    chain->traln->setFrequencySave(prevFrequRates[i], model, i); 

  chain->traln->initRevMat(model);
}

static void recordFrequRates(Chain *chain, int model, int numFrequRates, double *prevFrequRates)
{
  pInfo *partition = chain->traln->getPartition(model);

  assert(partition->dataType == DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    prevFrequRates[i] = chain->traln->getFrequency(model,i);
}



static void frequencySliderApply(Chain *chain, proposalFunction *pf)
{  
  int model = chain->getChainRand()->drawRandInt( chain->traln->getNumberOfPartitions()); 
  perPartitionInfo *info = pf->remembrance.partInfo; 
  pInfo*partition = chain->traln->getPartition(model); 

  int numFreq = partition->states; 
  info->modelNum = model; 
  info->numFreqs = numFreq; 

  recordFrequRates(chain, model, numFreq, info->frequencies); 
  
  int paramToChange = chain->getChainRand()->drawRandInt( numFreq); 
  double curv = chain->traln->getFrequency(model, paramToChange);   
  double newVal = fabs(chain->getChainRand()->drawFromSlidingWindow(curv, pf->parameters.slidWinSize)); 
  chain->traln->setFrequencySave(newVal,model, paramToChange); 

  double sum = 0;
  for(int i = 0; i < numFreq; ++i)
    sum += chain->traln->getFrequency(model,i); 
  for(int i = 0; i < numFreq; ++i)
    chain->traln->setFrequencySave(
				   chain->traln->getFrequency(model,i) / sum,
				   model,i); 

  
  chain->traln->initRevMat(model);
}



void frequency_proposal_apply(Chain * chain, proposalFunction *pf)
{
  int model  = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions());
  perPartitionInfo *info = pf->remembrance.partInfo; 

  pInfo *partition = chain->traln->getPartition( model); 

  int numFreq = partition->states; 
  info->modelNum = model;   
  info->numFreqs = numFreq; 

  recordFrequRates(chain, model, numFreq, info->frequencies);
  
  // double r[numFreq];  
  double* r  =  (double*)exa_calloc(numFreq, sizeof(double)); 
  for(int state = 0;state < numFreq ; state ++)
    {
      double curv = chain->traln->getFrequency(model,state);
      r[state] = chain->getChainRand()->drawRandBiUnif(curv);       
      chain->hastings*=curv/r[state];      
    }

  double sum=0;  
  for(int state = 0;state< numFreq ; state ++)    
    sum+=r[state]; 

  for(int state = 0; state< numFreq ; state ++)    
    chain->traln->setFrequencySave(r[state]/sum, model, state); 

  //recalculate eigens
  chain->traln->initRevMat(model);
  // exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
  exa_free(r); 

  /* chain->curprior=get_frequency_prior(chain, tr->partitionData[chain->frequRemem.model].frequencies); */
  /* chain->newprior=get_frequency_prior(chain, chain->frequRemem.curFrequRates); */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); */
}



void frequency_dirichlet_proposal_apply(Chain * chain, proposalFunction *pf)
{
   
  int model  = chain->getChainRand()->drawRandInt(chain->traln->getNumberOfPartitions());
  perPartitionInfo *info = pf->remembrance.partInfo; 

  pInfo *partition = chain->traln->getPartition( model); 

  int numFreq = partition->states; 
  info->modelNum = model;   
  info->numFreqs = numFreq; 

  recordFrequRates(chain, model, numFreq, info->frequencies);
  
  double* r  =  (double*)exa_calloc(numFreq, sizeof(double)); 
  
  
        
  //drawRandDirichlet(chain, r, partition->frequencies, 1.0, numFreq);
  chain->getChainRand()->drawDirichletExpected( r, partition->frequencies, pf->parameters.dirichletAlpha, numFreq);
  chain->hastings= densityDirichlet(partition->frequencies, r, numFreq) / densityDirichlet(r, partition->frequencies, numFreq);

  /*
  //Validation-----------------------------------------
  printf("frequencies: ");
  for(int state = 0;state < numFreq ; state ++)
  {
  printf("%f ", partition->frequencies[state]); 
  }
  printf("\n");
  printf("r:           ");
  for(int state = 0;state < numFreq ; state ++)
  {
  printf("%f ",r[state]); 
  }
  printf("\n");
  printf("hastings %f \n", chain->hastings);
  printf("\n");
  //--------------------------------------------------
  */
  
      
  for(int state = 0;state < numFreq ; state ++)
    {
      if (r[state]<FREQ_MIN)
	r[state]=FREQ_MIN+FREQ_MIN*FREQ_MIN*(numFreq-1);//with this minima sheme, minima should be satisfied even after rescaling.
    }
      
  double sum=0;  
  for(int state = 0;state< numFreq ; state ++)    
    sum+=r[state]; 

  for(int state = 0; state< numFreq ; state ++)    
    partition->frequencies[state]=r[state]/sum; 
  

  chain->traln->initRevMat(model);

  exa_free(r);
}

void frequency_proposal_reset(Chain * chain, proposalFunction *pf)
{
  perPartitionInfo *info = pf->remembrance.partInfo; 
  restore_frequ_rates(chain, info->modelNum, info->numFreqs, info->frequencies);
}



/*
 * should be sliding window proposal
 */

void edit_subs_rates(Chain *chain, int model, int subRatePos, double subRateValue)
{
  pInfo *partition = chain->traln->getPartition(model); 

  assert(partition->dataType == DNA_DATA);
  assert(subRateValue <= RATE_MAX && subRateValue >= RATE_MIN);
  int states = partition->states; 
  int numSubsRates = (states * states - states) / 2;
  assert(subRatePos >= 0 && subRatePos < numSubsRates);
  chain->traln->setSubstSave(subRateValue, model, subRatePos); 
}



static void simple_model_proposal_reset(Chain * chain, proposalFunction *pf)
{
  perPartitionInfo *info = pf->remembrance.partInfo; 
  restore_subs_rates(chain, info->modelNum, info->numRates, info->substRates);
}

static void branchLengthWindowApply(Chain *chain, proposalFunction *pf)
{  
  tree *tr = chain->traln->getTr(); 
  branch b =  chain->getChainRand()->drawBranchUniform(*(chain->traln)); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  pf->remembrance.modifiedPath->clear(); 
  pf->remembrance.modifiedPath->append( b);
  pf->remembrance.modifiedPath->at(0).length[0]  = chain->traln->getBranchLength( p,0); 

  double zOld = chain->traln->getBranchLength( p,0); 
  double win =  pf->parameters.slidWinSize ; 
  double realZ = branchLengthToReal(tr, zOld); 
  double blNew = branchLengthToInternal(tr,fabs(chain->getChainRand()->drawFromSlidingWindow( realZ, win))); 
  chain->traln->setBranchLengthSave(blNew, 0,p);
  // p->z[0] = p->back->z[0] = blNew; 
}

static void branchLengthMultiplierApply(Chain *chain, proposalFunction *pf)
{
  auto traln =  chain->traln; 
  tree *tr = chain->traln->getTr(); 
  branch b =  chain->getChainRand()->drawBranchUniform(*traln); 

  nodeptr p = findNodeFromBranch(tr, b); 

  pf->remembrance.modifiedPath->clear(); 
  pf->remembrance.modifiedPath->append( b);
  assert(pf->remembrance.modifiedPath->size() == 1 ); 
  pf->remembrance.modifiedPath->at(0).length[0]  = traln->getBranchLength( p,0); 

  double
    multiplier = chain->getChainRand()->drawMultiplier( pf->parameters.multiplier); 
  assert(multiplier > 0.); 
  
  /* TODO how do we do that wiht multiple bls per branch?  */
  assert(chain->traln->getNumBranches() == 1); 

  double oldZ = traln->getBranchLength( p,0); 
  double newZ = pow( oldZ, multiplier) ; 
  /* printf("\t\tBL_MULTI: %g * %g => %g\n", branchLengthToReal(tr, p->z[0]), multiplier, branchLengthToReal(tr, newZ));  */


#ifdef PRINT_MULT  
  cout  << setprecision(6) << "bl-multi: " << branchLengthToReal(tr,oldZ) <<   " * " << multiplier << " = "  << branchLengthToReal(tr, newZ) << endl; 
#endif


  /* according to lakner2008  */
  chain->hastings *= multiplier; 
 
  /* just doing it for one right here */
  traln->setBranchLengthSave(newZ, 0, p); 
}






/**
   @brief log tuning for a parameter 
   
   @return tuned parameter    
 */
double tuneParameter(int batch, double accRatio, double parameter, boolean inverse )
{  
  double delta = 1.0 / sqrt(batch);
  delta = 0.01 < delta ? 0.01 : delta;

  double logTuning = log(parameter);
  
  if(inverse)
    logTuning += (accRatio > TARGET_RATIO)  ? -delta : +delta ;
  else 
    logTuning += (accRatio > TARGET_RATIO)  ? +delta : -delta ;
  
  double newTuning = exp(logTuning);

  /* TODO min+max tuning?  */
  
  double minTuning = 1e-8,
    maxTuning = 1e5; 
  if (minTuning <  newTuning && newTuning < maxTuning)
    return  newTuning; 
  else 
    return parameter; 
}






/**
   @brief tunes the BL multiplier
 */ 
static void autotuneMultiplier( proposalFunction *pf, SuccessCtr *ctr)
{
  double *parameter = &(pf->parameters.multiplier); 

  double newParam = tuneParameter(ctr->getBatch(), ctr->getRatioInLastInterval(), *parameter, FALSE); 

#ifdef DEBUG_PRINT_TUNE_INFO
  cout << pf->name << ": with ratio " << ctr->getRatioInLastInterval() << ": "<< ((newParam < *parameter ) ? "reducing" : "increasing") <<  "\t" << *parameter << "," << newParam << endl; 
#endif 

  *parameter = newParam; 
  ctr->nextBatch();
}


/**
   @brief autotunes sliding windows 
   
*/ 
static void autotuneSlidingWindow(proposalFunction *pf, SuccessCtr *ctr)
{
  double *parameter = &(pf->parameters.slidWinSize); 
  double newParam = tuneParameter(ctr->getBatch() , ctr->getRatioInLastInterval(), *parameter, FALSE  ); 
  
#ifdef DEBUG_PRINT_TUNE_INFO
  cout << pf->name << ": with ratio " << ctr->getRatioInLastInterval() << ": "<< ((newParam < *parameter ) ? "reducing" : "increasing") <<  "\t" << *parameter << "," << newParam << endl; 
#endif

  *parameter = newParam; 
  ctr->nextBatch();
}




/** @brief tune the alpha of a dirichlet distribution */ 
static void autotuneDirichletAlpha(proposalFunction *pf, SuccessCtr *ctr)
{
  double *parameter = &(pf->parameters.dirichletAlpha); 

  double newParam = tuneParameter(ctr->getBatch(), ctr->getRatioInLastInterval(), *parameter, TRUE); 

#ifdef DEBUG_PRINT_TUNE_INFO
  cout << pf->name << ": with ratio " << ctr->getRatioInLastInterval() << ": "<< ((newParam < *parameter ) ? "reducing" : "increasing") <<  "\t" << *parameter << "," << newParam << endl; 
#endif

  * parameter = newParam; 
  ctr->nextBatch();
}


void branchLengthReset(Chain *chain, proposalFunction *pf)
{
  assert(pf->remembrance.modifiedPath->size() == 1 ); 
  branch b = pf->remembrance.modifiedPath->at(0); 
  pf->remembrance.modifiedPath->clear(); 
  tree *tr = chain->traln->getTr(); 
  nodeptr p = findNodeFromBranch(tr, b); 
  chain->traln->setBranchLengthSave(b.length[0], 0,p); 
  // p->z[0] = p->back->z[0] = b.length[0]; 
}


void guided_branch_length_proposal_apply(Chain *chain, proposalFunction *pf)
{
  // DUMMY, just for making it compile. 

  // please delete this function, once you have implemented it or added the source
}


void initProposalFunction( proposal_type type, vector<double> weights, proposalFunction **result)
{
  if(weights[type] == 0)
    {
      *result = NULL; 
      return ; 
    }  

  *result = (proposalFunction*)exa_calloc(1,sizeof(proposalFunction));   

  proposalFunction *ptr = *result; 
  ptr->ptype = (proposal_type)type; 
  ptr->initWeight = weights[type]; 
  ptr->currentWeight = ptr->initWeight; 
  ptr->relativeWeight = weights[type] ; 


  /* 
     TODO@kassian
     some proposals may be broken      
  */

  switch(type)
    {
    case UPDATE_SINGLE_BL: 
      ptr->eval_lnl = evalBranch;
      ptr->autotune = autotuneSlidingWindow; 
      ptr->remembrance.modifiedPath = new Path(); 
      // createStack(&(ptr->remembrance.modifiedPath)); 
      ptr->apply_func = branchLengthWindowApply; 
      ptr->reset_func =  branchLengthReset;
      ptr->category = BRANCH_LENGTHS; 
      ptr->parameters.slidWinSize = INIT_BL_SLID_WIN; 
      ptr->name = "singleBLSlidWin"; 
      break; 
    case UPDATE_SINGLE_BL_EXP: 
      ptr->eval_lnl = dummy_eval;
      ptr->remembrance.topoRec = (topoRecord*)exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func	=  random_branch_length_proposal_apply;
      ptr->reset_func =  random_branch_length_proposal_reset;
      ptr->category = BRANCH_LENGTHS; 
      ptr->name  = "singleBlExp"; 
      break;
    case UPDATE_SINGLE_BL_GUIDED: 
      ptr->eval_lnl = evalBranch;
      ptr->remembrance.modifiedPath = new Path(); 
      // createStack(&(ptr->remembrance.modifiedPath));       
      ptr->apply_func =  guided_branch_length_proposal_apply;
      ptr->reset_func =  branchLengthReset;
      ptr->category = BRANCH_LENGTHS; 
      ptr->name  = "singleBlGuided"; 
      break; 
    case UPDATE_SINGLE_BL_BIUNIF: 
      ptr->eval_lnl = dummy_eval;
      ptr->remembrance.topoRec = (topoRecord*)exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func	=  random_branch_length_proposal_apply;
      ptr->reset_func =  random_branch_length_proposal_reset;
      ptr->name = "singleBLBiunif"; 
      ptr->category = BRANCH_LENGTHS; 
      break; 
    case UPDATE_MODEL_SINGLE_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func	=  single_biunif_model_proposal_apply;
      ptr->reset_func =  simple_model_proposal_reset;
      ptr->name = "singleModelBiunif"; 
      ptr->category = SUBSTITUTION_RATES; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      break; 
    case UPDATE_MODEL_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->name = "modelBiunif"; 
      ptr->category = SUBSTITUTION_RATES; 
      ptr->apply_func = single_biunif_model_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      break; 
    case UPDATE_MODEL_ALL_BIUNIF:
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = all_biunif_model_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->name = "modelAllBiunif"; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->category = SUBSTITUTION_RATES; 
      break; 
    case UPDATE_FREQUENCIES_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = frequency_proposal_apply; 
      ptr->reset_func = frequency_proposal_reset; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "freqBiunif"; 
      ptr->category = FREQUENCIES; 
      break;
    case UPDATE_MODEL_PERM_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = NULL; 	/* TODO */
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "modelPermBiunif"; 
      ptr->category = SUBSTITUTION_RATES; 
      break;
    case BRANCH_LENGTHS_MULTIPLIER: 
      ptr->eval_lnl = evalBranch;
      ptr->autotune = autotuneMultiplier; 
      ptr->apply_func = branchLengthMultiplierApply; 
      ptr->reset_func = branchLengthReset; 
      ptr->remembrance.modifiedPath =  new Path(); 
      // createStack(&(ptr->remembrance.modifiedPath));
      ptr->parameters.multiplier = INIT_BL_MULT; 
      ptr->name = "branchMult"; 
      ptr->category = BRANCH_LENGTHS; 
      break; 
      /* TODO re-install PROPOSALADD anchor for script   */
    default : 
      {
	printf("unknown value %d\n",type ); 
	assert(0) ; 
      }
    }  
}
