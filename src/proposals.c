/**
   @file proposals.c

   @brief All proposals that make the chains of ExaBayes move in its
   space.
    
 */ 



#include "common.h"

#include "axml.h"
#include "proposalStructs.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"
#include "convergence.h"
#include "chain.h"
#include "topology-utils.h" 
#include "eval.h"
#include "adapters.h"
#include "exa-topology.h"


void expensiveVerify(tree *tr); 

nodeptr select_random_subtree(state *chain, tree *tr);
void edit_subs_rates(state *chain, int model, int subRatePos, double subRateValue);



//reads proposalWeights p and sets t for logistic function such that f(t)=p
void findLogisticT(state *chain){
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
void findLogisticP(state *chain){
  static double leftShift=2.0;
  for(int i = 0; i < NUM_PROPOSALS;++i){
    if(chain->proposalLogisticT[i]>=0){
    chain->proposalWeights[i]=(1.0/(1.0+exp(leftShift-chain->proposalLogisticT[i])) );  
  //  printf("p(i): %f e^..: %f t: %f\n",chain->proposalWeights[i], (1.0/(1.0+exp(leftShift-chain->proposalLogisticT[i])) ), chain->proposalLogisticT[i]);
   }
  }
}



static void recordSubsRates(state *chain, int model, int numSubsRates, double *prevSubsRates)
{
  pInfo *partition = getPartition(chain,model); 
  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = partition->substRates[i];
}


static void reset_branch_length(nodeptr p, int numBranches)
{
  int i;
  for(i = 0; i < numBranches; i++)
    {
      assert(p->z_tmp[i] == p->back->z_tmp[i]);
      p->z[i] = p->back->z[i] = p->z_tmp[i];   /* restore saved value */
    }
}


static void record_branch_info(nodeptr p, double *bl, int numBranches)
{
  int i;
  for(i = 0; i < numBranches; i++)
    bl[i] = p->z[i];
}


/* static node *extended_spr_traverse( tree *tr, node *n ) { */

/*   double randprop = drawRandDouble(); */
  
/*   if( isTip(n->number, tr->mxtips ) || randprop < 0.5 ) { */
/*     return n; */
/*   } else if( randprop < 0.75 ) { */
/*     spr_depth++; */
/*     return extended_spr_traverse( tr, n->next->back ); */
/*   } else { */
/*     spr_depth++; */
/*     return extended_spr_traverse( tr, n->next->next->back ); */
/*   } */

/* } */




/* static node *extended_spr( tree *tr, node *n )  */
/* { */
/*   if( isTip(n->number, tr->mxtips ) )  */
/*     return n; */

/*   double randprop = drawRandDouble(); */
/*   if( randprop < 0.5 )  */
/*     return extended_spr_traverse( tr, n->next->back ); */
/*   else  */
/*     return extended_spr_traverse( tr, n->next->next->back ); */
/* } */


int radius = 0; 

int extended_spr_traverse(state *chain, nodeptr *insertNode)
{
  int 
    result, 
    r = drawRandDouble01(chain); 
  
  radius++;   

  /* *insertNode = NULL;  */
  if(r < 0.5 )
    {


       *insertNode = isTip((*insertNode)->next->number, chain->tr->mxtips) ? (*insertNode)->next :  (*insertNode)->next->back; 

     // *insertNode = (*insertNode)->next->back ;
    }
  else 
    {
       *insertNode = isTip((*insertNode)->next->next->number, chain->tr->mxtips) ? (*insertNode)->next->next :  (*insertNode)->next->next->back; 
      //*insertNode = (*insertNode)->next->next->back;
    }

  /*
  if(isTip( (*insertNode)->number, chain->tr->mxtips))
  {
     if(randNum == 0)
    direction--;
    else
    direction++;
  }*/
  /*
  if(direction % 3 == 0)
    {
      if(randNum == 0)
	*insertNode = (*insertNode)->next->back; 
      else 
	*insertNode = (*insertNode)->next->next->back; 
    }
  else if(direction % 3 == 1)
    {
      if(randNum == 0)
	*insertNode = (*insertNode)->next->next->back; 	
      else 
	*insertNode = (*insertNode)->back; 
    }
  else 
    {
      if(randNum == 0)
	*insertNode = (*insertNode)->back; 
      else 
	*insertNode = (*insertNode)->next->back; 
    }
    
     if(isTip( (*insertNode)->number, chain->tr->mxtips))
  {
     if(randNum == 0)
    direction--;
    else
    direction++;
  }
  */
  double randprop = drawRandDouble01(chain);
  result = randprop < chain->eSprStopProb; 
  
  return result; 
}








static void extended_spr_apply(state *chain, int pSubType)
{
  tree * tr = chain->tr;

#ifdef DEBUG_SHOW_TREE
  char tmp[10000]; 
  Tree2stringNexus(tmp, tr, tr->start->back, 0); 
  if(processID==0)
    printf("topo before: %s\n", tmp);
#endif

  
  nodeptr    
    p = select_random_subtree(chain,tr);

#if 0
  parsimonySPR(p, tr);
#endif  

  chain->sprMoveRemem.p = p;
  chain->sprMoveRemem.nb  = p->next->back;
  chain->sprMoveRemem.nnb = p->next->next->back;
  
  record_branch_info(chain->sprMoveRemem.nb, chain->sprMoveRemem.nbz, getNumBranches(chain->tr));
  record_branch_info(chain->sprMoveRemem.nnb, chain->sprMoveRemem.nnbz, getNumBranches(chain->tr));
  

  /* initial remapping of BL of nodes adjacent to pruned node  */
  double zqr[NUM_BRANCHES];
  
  chain->hastings = 1; 
  
  for(int i = 0; i < getNumBranches(tr); i++)
    {
      if(pSubType == STANDARD)
	{
	  zqr[i] = chain->sprMoveRemem.nb->z[i] * chain->sprMoveRemem.nnb->z[i];
	
	  chain->hastings *= log(zqr[i]);
	}
      else if(pSubType == SPR_MAPPED)
	{
	  zqr[i] = chain->sprMoveRemem.nb->z[i] ; 
	}
      else if(pSubType == SPR_ADJUST)
	{
	  zqr[i] = chain->sprMoveRemem.nb->z[i] * chain->sprMoveRemem.nnb->z[i];  
	  chain->hastings *= log(zqr[i]);
	  zqr[i] = sqrt(zqr[i]);
	}
      if(zqr[i] > zmax) zqr[i] = zmax;
      if(zqr[i] < zmin) zqr[i] = zmin;
    }

  hookup(chain->sprMoveRemem.nb, chain->sprMoveRemem.nnb, zqr, getNumBranches(tr)); 
  p->next->next->back = p->next->back = (node *) NULL;
  /* done remove node p (omitted BL opt) */

  nodeptr initNode = NULL; 
  nodeptr curNode = NULL; 
  boolean remapBL = FALSE; 
  if(drawRandDouble01(chain) > 0.5)
    {
      curNode = chain->sprMoveRemem.nb; 
      remapBL = TRUE ; 
    }
  else 
    {
      curNode = chain->sprMoveRemem.nnb; 
    }
  initNode = curNode; 

  int accepted = FALSE;   

  radius = 0; 
//printf("curNode:  back next nextnext\n");
  while( NOT  accepted)
    {
      accepted = extended_spr_traverse(chain, &curNode); 
  //    if(processID == 0) 
//	printf("%d:  %d  %d  %d\n",curNode->number, curNode->back->number, curNode->next->back->number, curNode->next->next->back->number); 

      /* needed for spr remap */
      if(curNode == initNode)
	remapBL = NOT remapBL; 
    }



  chain->newprior = 1; 
  chain->curprior = 1; 

  chain->sprMoveRemem.q = curNode;
  chain->sprMoveRemem.r = chain->sprMoveRemem.q->back;
  record_branch_info(chain->sprMoveRemem.q, chain->brLenRemem.qz, getNumBranches(chain->tr));

  
  
   if(pSubType == STANDARD){
     for(int branchCount=0; branchCount< getNumBranches(chain->tr); branchCount++)
     {
      chain->hastings/=log(curNode->z[branchCount]);
     }
   }else if(pSubType == SPR_ADJUST){
     for(int branchCount=0; branchCount< getNumBranches(chain->tr); branchCount++) /*  */
     {
      chain->hastings/=(2*log(curNode->z[branchCount]));
     }
   }
  
 // Thorough = 0;
  // assert(tr->thoroughInsertion == 0);
  /* insertBIG wont change the BL if we are not in thorough mode */

  if(pSubType == STANDARD)
    {
      insertWithUnifBL(chain->sprMoveRemem.p, chain->sprMoveRemem.q, getNumBranches(chain->tr));
  /*   for(int branchCount=0; branchCount<chain->tr->numBranches; branchCount++)
     {
       chain->hastings/=log(chain->brLenRemem.qz[branchCount]);
     }
     */
   //  printf("hastings: %f\n", chain->hastings);
      // insertBIG(chain->tr, chain->sprMoveRemem.p, chain->sprMoveRemem.q, chain->tr->numBranches);
    }
  else if(pSubType == SPR_MAPPED)
    {
      double *neighborZ = remapBL ? chain->sprMoveRemem.nbz :  chain->sprMoveRemem.nnbz; 
      
        if( remapBL ) 
	  {
	    for(int i = 0; i < getNumBranches(tr); ++i)
	      chain->sprMoveRemem.nb->z[i] = chain->sprMoveRemem.nb->back->z[i] = chain->sprMoveRemem.nnbz[i]; 
	  }
      
	insertWithGenericBL(chain->sprMoveRemem.p, chain->sprMoveRemem.q, chain->sprMoveRemem.p->z, curNode->z, neighborZ, getNumBranches(tr));


      /* IMPORTANT TODO verify, that the mapping actually works, as we
	 had this in mind. Use topo print functions for that.
       */


      /* TODO PERFORMANCE */
      /* newviewGeneric(tr, chain->sprMoveRemem.p, FALSE); */
    }
     else if(pSubType == SPR_ADJUST)
    {
      insertWithUnifBLScaled(chain->sprMoveRemem.p, chain->sprMoveRemem.q, 2.0,  getNumBranches(chain->tr));
    }

//   if( (pSubType == SPR_MAPPED ) && remapBL) 
//     {
//       for(int i = 0; i < tr->numBranches; ++i)
// 	chain->sprMoveRemem.nb->z[i] = chain->sprMoveRemem.nb->back->z[i] = chain->sprMoveRemem.nnbz[i]; 
//     }

  /* TODO problem here? is not tip?? FIXED: should never be a tip anymore. Possible since we are interested in edges and each edge connects to at least one inner node*/
#if 0 
  evaluateGeneric(chain->tr, chain->sprMoveRemem.p->next->next, FALSE);
#else   
  evaluateGenericWrapper(chain, tr->start, TRUE);
#endif
}





static void extended_spr_reset(state * chain)
{
  tree *tr = chain->tr; 

  /* prune the insertion */
  hookup(chain->sprMoveRemem.q, chain->sprMoveRemem.r, chain->brLenRemem.qz, getNumBranches(chain->tr));

  chain->sprMoveRemem.p->next->next->back = chain->sprMoveRemem.p->next->back = (nodeptr) NULL;
  /*  */
  /* insert the pruned tree in its original node */
  hookup(chain->sprMoveRemem.p->next,        chain->sprMoveRemem.nb, chain->sprMoveRemem.nbz, getNumBranches(chain->tr));
  hookup(chain->sprMoveRemem.p->next->next, chain->sprMoveRemem.nnb, chain->sprMoveRemem.nnbz, getNumBranches(chain->tr));
  
  if(processID == 0)
    {

#ifdef DEBUG_SHOW_TREE
      char tmp[100000];
      Tree2stringNexus(tmp, chain->tr, chain->tr->start->back, 0); 
      printf("topo reset: %s\n", tmp); 
#endif
    }

  evaluateGenericWrapper(chain, tr->start, TRUE);
  
  exa_newViewGeneric(chain, chain->sprMoveRemem.p, FALSE); 
  double val1 = chain->tr->likelihood; 
  
  exa_newViewGeneric(chain, chain->tr->start, TRUE);
  double  val2 = chain->tr->likelihood; 

  assert( fabs ( val2 - val1 ) < 0.0001 ); 

}


//--------Alpha-Proposal-for-GAMMA-----------------------------------------------------------------

double get_alpha_prior(state *chain )
{
  
 return 1;//TODO obviously needs acctual prior 
}





static void simple_gamma_proposal_apply(state * chain, int pSubType)
{
  tree *tr = chain->tr; 

  assert(tr != NULL); 

  //TODO: add safety to max and min values
  double newalpha, curv, r,mx,mn;
  chain->modelRemem.model=drawRandInt(chain, getNumberOfPartitions(tr));

  pInfo *partition = getPartition(chain,chain->modelRemem.model);

  curv = partition->alpha;
  chain->gammaRemem.curAlpha = curv;


  switch(pSubType)
        {
        case STANDARD://simple sliding window
	  r = drawRandDouble01(chain);
	  mn = curv-(chain->gammaRemem.gm_sliding_window_w/2);
	  mx = curv+(chain->gammaRemem.gm_sliding_window_w/2);
	  newalpha = fabs(mn + r * (mx-mn));
	  break;
        case EXP_DISTR:
          newalpha  = drawRandExp(chain,1/curv);
          break;
        default:
          assert(0);
        }
  /* Ensure always you stay within this range */
  if(newalpha > ALPHA_MAX) newalpha = ALPHA_MAX;
  if(newalpha < ALPHA_MIN) newalpha = ALPHA_MIN;
  
  chain->hastings = 1; //since it is symmetrical, hastings=1
  chain->newprior = get_alpha_prior(chain); 
  chain->curprior = get_alpha_prior(chain); 
  
  partition->alpha = newalpha;

  makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);

  evaluateOnePartition(chain, tr->start, TRUE, chain->modelRemem.model); 
}


// static void exp_gamma_proposal_apply(state * chain, int pSubType)
// {
//   double newalpha, curv;
//   curv = chain->tr->partitionData[chain->modelRemem.model].alpha;
//   chain->gammaRemem.curAlpha = curv;
//   newalpha  = drawRandExp(1/curv);
//   
// 
//  
//   /* Ensure always you stay within this range */
//   while(newalpha > ALPHA_MAX || newalpha < ALPHA_MIN){
//   if(newalpha > ALPHA_MAX) newalpha = 2*ALPHA_MAX-newalpha;
//   if(newalpha < ALPHA_MIN) newalpha = 2*ALPHA_MIN-newalpha;
//   }
//   chain->hastings = (1/newalpha)*exp(-(1/newalpha)*curv)/((1/curv)*exp(-(1/curv)*newalpha)); //TODO do not ignore reflection
//   /* chain->newprior = get_alpha_prior(chain);  */
//   /* chain->curprior = get_alpha_prior(chain);  */
//   
//   chain->tr->partitionData[chain->modelRemem.model].alpha = newalpha;
//   
//   makeGammaCats(chain->tr->partitionData[chain->modelRemem.model].alpha, chain->tr->partitionData[chain->modelRemem.model].gammaRates, 4, chain->tr->useMedian);
// 
//   evaluateGeneric(chain->tr, chain->tr->start, TRUE);
// }


static void simple_gamma_proposal_reset(state * chain)
{
  tree *tr = chain->tr; 
  pInfo *partition = getPartition(chain, chain->modelRemem.model) ; 
  
  partition->alpha = chain->gammaRemem.curAlpha; 

  makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);

  evaluateOnePartition(chain, tr->start, TRUE, chain->modelRemem.model); 
}

//------------------------------------------------------------------------------

/* TODO we do not even use this function, do we? NOTE Now we do ;) */
void penalize(state *chain, int which_proposal, int acceptance)
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





//NOTE: should only be called at the very beginning. Afterwards a probability sum of 1.0 is not required.

void normalizeProposalWeights(state *chain)
{
  double sum = 0 ; 
  for(int i = 0; i < NUM_PROPOSALS;++i)
    sum += chain->proposalWeights[i]; 
  
  for(int i = 0; i < NUM_PROPOSALS;++i)
    chain->proposalWeights[i] /= sum; 
  
  findLogisticT(chain);
}



static void simple_model_proposal_apply(state *chain, int pSubType)//llpqr
{
  tree *tr = chain->tr; 
  
  //TODO: add safety to max and min values
  //record the old ones
  chain->modelRemem.model=drawRandInt(chain, getNumberOfPartitions(chain->tr));
  recordSubsRates(chain, chain->modelRemem.model, chain->modelRemem.numSubsRates, chain->modelRemem.curSubsRates);

  pInfo *partition = getPartition(chain, chain->modelRemem.model); 

  //choose a random set of model params,
  //probably with dirichlet proposal
  //with uniform probabilities, no need to have other
  int state, changeState = 0;
  double new_value = 0,curv;
  double r;

  
  double mx,mn; //for sliding window
  
  int list[chain->modelRemem.numSubsRates];//for biunif_distr and biunif_perm_distr
  
  /* int numberOfEdits;//for biunif_perm_distr */
  
  switch(pSubType)
        {
	  case STANDARD:
	    case BIUNIF_DISTR:
	      /* numberOfEdits=chain->modelRemem.numSubsRates; */
	      break;
	      case BIUNIF_PERM_DISTR:
		/* numberOfEdits=drawRandInt(chain, chain->modelRemem.numSubsRates); */
		break;
	} 
  
  
  chain->hastings=1.0;
  
    for(state = 0;state<chain->modelRemem.numSubsRates ; state ++)
    {
      
      switch(pSubType)
        {
        case STANDARD: //using the branch length sliding window for a test    
	  changeState=state;
	  curv = partition->substRates[state];
	  r =  drawRandDouble01(chain);
	  mn = curv-(chain->modelRemem.rt_sliding_window_w/2);
	  mx = curv+(chain->modelRemem.rt_sliding_window_w/2);
	
	  new_value = fabs(mn + r * (mx-mn));
	  
	  /* Ensure always you stay within this range */
      if(new_value > RATE_MAX) new_value = RATE_MAX;
      if(new_value < RATE_MIN) new_value = RATE_MIN;
      
	  break;
        
	case BIUNIF_DISTR:
	  changeState=drawRandInt(chain, chain->modelRemem.numSubsRates);
      if(list[changeState]!=1)
      {
      list[changeState]=1;;      
      curv = partition->substRates[changeState];
      r =  drawRandBiUnif(chain, curv);
      new_value = r;
      chain->hastings*=curv/new_value;
 
      while(new_value> RATE_MAX|| new_value< RATE_MIN)
      {
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      }
      break;
      case BIUNIF_PERM_DISTR://TODO NOT used. Lower values than before (without subType)
	drawPermutation(chain, list, chain->modelRemem.numSubsRates);
	 curv = partition->substRates[list[state]];
	 r =  drawRandBiUnif(chain, curv);
       new_value = r;
      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       chain->hastings*=curv/new_value;
	 break;
	 
	case SINGLE_BIUNIF://TODO not used. Figure out error
	  drawPermutation(chain,list, chain->modelRemem.numSubsRates);
	 curv = partition->substRates[list[state]];
	 r =  drawRandBiUnif(chain, curv);
       new_value = r;
      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       chain->hastings*=curv/new_value;
	 break;
	 
        default:
          assert(0);
        
      

    }//end of switch
      edit_subs_rates(chain,chain->modelRemem.model, changeState, new_value);
    }
  //recalculate eigens

    exa_initReversibleGTR(chain, chain->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
    
  /* TODO: need to broadcast rates here for parallel version ! */

  evaluateOnePartition(chain, tr->start, TRUE, chain->modelRemem.model); /* 2. re-traverse the full tree to update all vectors */

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


static void perm_biunif_model_proposal_apply(state *chain, int pSubType)
{
  tree *tr = chain->tr; 

  //record the old one 
  chain->modelRemem.model=drawRandInt(chain,getNumberOfPartitions(chain->tr));
  
  pInfo *partition = getPartition(chain , chain->modelRemem.model) ; 

  recordSubsRates(chain, chain->modelRemem.model, chain->modelRemem.numSubsRates, chain->modelRemem.curSubsRates);
  int state, randNumber;
  double new_value,curv;
  double r;
  
  randNumber=drawRandInt(chain,chain->modelRemem.numSubsRates);
    int perm[chain->modelRemem.numSubsRates];
    drawPermutation(chain,perm, chain->modelRemem.numSubsRates);
  
   chain->hastings=1.0;
  for(state = 0;state<randNumber ; state ++)
    {           
      curv = partition->substRates[perm[state]];
      r =  drawRandBiUnif(chain,curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       chain->hastings*=curv/new_value;
      edit_subs_rates(chain,chain->modelRemem.model, perm[state], new_value);
      
    }
      
  

  exa_initReversibleGTR(chain, chain->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  evaluateOnePartition(chain, tr->start, TRUE, chain->modelRemem.model); /* 2. re-traverse the full tree to update all vectors */  
}


static void single_biunif_model_proposal_apply(state *chain,int pSubType)//NOTE whenever a model parameter changes, all branch lengths have to be re-normalized with 1/fracchange. Additionally we always must do a full tree traversal to get the likelihood. So updating a single parameter is rather expensive, .
{
  tree *tr = chain->tr; 
  //record the old one //TODO sufficient to store single value.
  chain->modelRemem.model=drawRandInt(chain,getNumberOfPartitions(chain->tr)); 
  
  pInfo *partition = getPartition(chain,chain->modelRemem.model) ; 

  recordSubsRates(chain, chain->modelRemem.model, chain->modelRemem.numSubsRates, chain->modelRemem.curSubsRates);
  //choose a random set parameter,
  //with uniform probabilities

  int  randState=drawRandInt(chain,chain->modelRemem.numSubsRates);

  double new_value,curv;
  double r;
  
  //int state=drawRandInt(chain->modelRemem.numSubsRates);
  

  curv = partition->substRates[randState];
  r =  drawRandBiUnif(chain,curv);

  new_value = r;
      
  while(new_value> RATE_MAX|| new_value< RATE_MIN){
    if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
    if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
  }

  edit_subs_rates(chain,chain->modelRemem.model, randState, new_value);

  chain->hastings=curv/new_value;

  exa_initReversibleGTR(chain, chain->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
  
  evaluateOnePartition(chain, tr->start, TRUE, chain->modelRemem.model); /* 2. re-traverse the full tree to update all vectors */
}

static void all_biunif_model_proposal_apply(state *chain, int pSubType)
{
  tree *tr = chain->tr; 
  
  //record the old one 
  chain->modelRemem.model=drawRandInt(chain,getNumberOfPartitions(chain->tr));
  pInfo *partition = getPartition(chain, chain->modelRemem.model) ; 

  recordSubsRates(chain, chain->modelRemem.model, chain->modelRemem.numSubsRates, chain->modelRemem.curSubsRates);
  //choose a random set parameter,
  //with uniform probabilities
  int state;
  double new_value,curv;
  double r;
  
  
   chain->hastings=1.0;

  for(state = 0;state<chain->modelRemem.numSubsRates ; state ++)
    {
      curv = partition->substRates[state]; 
      r =  drawRandBiUnif(chain,curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       chain->hastings*=curv/new_value;
      edit_subs_rates(chain,chain->modelRemem.model, state, new_value);
    }

  exa_initReversibleGTR(chain, chain->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  evaluateOnePartition(chain, tr->start, TRUE, chain->modelRemem.model); /* 2. re-traverse the full tree to update all vectors */
}

static void restore_subs_rates(state *chain, analdef *adef, int model, int numSubsRates, double *prevSubsRates)
{
  tree *tr = chain->tr; 

  pInfo *partition = getPartition(chain, model); 

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)	
    partition->substRates[i] = prevSubsRates[i]; 

  exa_initReversibleGTR(chain, model);

  /* TODO need to broadcast rates here for parallel version */

  evaluateOnePartition(chain, tr->start, TRUE, model); 
}


//--------Branch-Length_Proposals---------------------------------------------------

double get_branch_length_prior( state *chain)
{//TODO decide on sensible prior
  return 1;  
}

//setting this out to allow for other types of setting
static void set_branch_length_sliding_window(state *chain, nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
{
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
	  p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /* keep current value */
	}
      r = drawRandDouble01(chain);

      real_z = -log(p->z[i]) * s->tr->fracchange;
      
      //     printf( "z: %f %f\n", p->z[i], real_z );
    
      mn = real_z-(s->brLenRemem.bl_sliding_window_w/2);
      mx = real_z+(s->brLenRemem.bl_sliding_window_w/2);
      newZValue = exp(-(fabs(mn + r * (mx-mn)/s->tr->fracchange )));
    
      
      s->hastings=1;
      
      /* newValue = mn + r * (mx-mn); */
      /* s->newprior=get_branch_length_prior(&newValue); */
    
      /* Ensure always you stay within this range */
      if(newZValue > zmax) newZValue = zmax;
      if(newZValue < zmin) newZValue = zmin;
      assert(newZValue <= zmax && newZValue >= zmin);
      //     printf( "z: %f %f %f\n", p->z[i], newZValue, real_z );
      p->z[i] = p->back->z[i] = newZValue;
      //assuming this will be visiting each node, and multiple threads won't be accessing this
      //s->bl_prior += log(exp_pdf(s->bl_prior_exp_lambda,newZValue));
      //s->bl_prior += 1;
    }
}

static void set_branch_length_biunif(state *chain, nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
{
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
	  p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /* keep current value */
	}
      
      real_z = -log(p->z[i]) * s->tr->fracchange; //convert from exponential to real form
      r = drawRandBiUnif(chain, real_z);//draw from [real_z/2,2*real_z]

    
      newZValue = exp(-(fabs( r/s->tr->fracchange )));//convert to exponential form
      
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
      p->z[i] = p->back->z[i] = newZValue;

    }
}

static void set_branch_length_exp(state *chain, nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
{
  
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
	  p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /* keep current value */
	}
   //   r = drawRandExp(lambda);
       real_z = -log(p->z[i]) * s->tr->fracchange;
       r = drawRandExp(chain, 1.0/real_z);
  
	s->hastings=(1/r)*exp(-(1/r)*real_z)/((1/real_z)*exp(-(1/real_z)*r));

	/* s->newprior=get_branch_length_prior(&r); */
	
      newZValue = exp(-(fabs(r /s->tr->fracchange )));
   
    
      /* Ensure always you stay within this range, reflect if necessary *///TODO reflects transformed bl, needs to reflect actual bl
      while(newZValue > zmax || newZValue < zmin){
      if(newZValue > zmax) newZValue = 2*zmax-newZValue;
      if(newZValue < zmin) newZValue = 2*zmin-newZValue;
      }
      assert(newZValue <= zmax && newZValue >= zmin);
      p->z[i] = p->back->z[i] = newZValue;
     }
     
}


static node *select_branch_by_id_dfs_rec( node *p, int *cur_count, int target, state *s ) {
  if( (*cur_count) == target ) 
    {
      return p;
    }
  
  if( !isTip( p->number, s->tr->mxtips )) {
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


static node *select_branch_by_id_dfs( node *p, int target, state *s ) {
  const int num_branches = (2 * s->tr->mxtips) - 3;
  
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
    
    if( isTip( q->number, s->tr->mxtips ) ) {
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




static void random_branch_length_proposal_apply(state * chain, int pSubType)
{
  int numBranches = getNumBranches(chain->tr);
  const int num_branches = (chain->tr->mxtips * 2) - 3;
  int target_branch = drawRandInt(chain,num_branches); 
  node *p = select_branch_by_id_dfs( chain->tr->start, target_branch, chain );
  
    switch(pSubType)
        {
        case STANDARD://simple sliding window
	     //for each branch get the current branch length
	    //pull a uniform like
	    //x = current, w =window
	    //uniform(x-w/2,x+w/2)
	  set_branch_length_sliding_window(chain,p, numBranches, chain, TRUE);
	  break;
        case EXP_DISTR:
	  set_branch_length_exp(chain,p, numBranches, chain, TRUE);
          break;
	case BIUNIF_DISTR:
	  set_branch_length_biunif(chain, p, numBranches, chain, TRUE);
	  break;
        default:
          assert(0);
        }


  chain->brLenRemem.single_bl_branch = target_branch;
  evaluateGenericWrapper(chain, p, FALSE); 
  //evaluateGeneric(chain->tr, chain->tr->start, TRUE); /* update the tr->likelihood *//FALSE seems to work
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

static void random_branch_length_proposal_reset(state * chain)
{
  node *p;
  assert( chain->brLenRemem.single_bl_branch != -1 );
  
  // ok, maybe it would be smarter to store the node ptr for rollback rather than re-search it...
  p = select_branch_by_id_dfs( chain->tr->start, chain->brLenRemem.single_bl_branch, chain );
  
  reset_branch_length(p, getNumBranches(chain->tr));
  //   printf( "reset bl: %p %f\n", p, p->z[0] );
  //update_all_branches(chain, TRUE);

  /* i am not so sure, if we really can do the FALSE here ; will raxml recognize that a branch length has changed? disabling it for now... 

     TODO I think, we should evaluate at the respctive node 
   */
#if 0 
  evaluateGenericWrapper(chain, chain->tr->start, FALSE);
#else 
  evaluateGenericWrapper(chain, chain->tr->start, TRUE );
#endif

 // evaluateGeneric(chain->tr, p, FALSE); //This yields a very slight likelihood difference.NOTE if we want exact likelihoods as before the proposal, we must evaluate from chain->tr->start, that is: evaluateGeneric(chain->tr, chain->tr->start, TRUE);
  chain->brLenRemem.single_bl_branch = -1;
}

double get_frequency_prior(state * chain)
{
 return 1; 
}

static void restore_frequ_rates(state *chain, analdef *adef, int model, int numFrequRates, double *prevFrequRates)
{
  tree *tr = chain->tr; 

  /* NOTICE: this function should not be called repeatedly  */

  pInfo *partition = getPartition(chain,model);

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    partition->frequencies[i] = prevFrequRates[i];

  exa_initReversibleGTR(chain, model);

  evaluateOnePartition(chain, tr->start, TRUE, model); 
}

static void recordFrequRates(state *chain, int model, int numFrequRates, double *prevFrequRates)
{
  pInfo *partition = getPartition(chain,model); 

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    prevFrequRates[i] = partition->frequencies[i];
}

void frequency_proposal_apply(state * chain, int pSubType)
{
  tree *tr = chain->tr; 

  chain->frequRemem.model=drawRandInt(chain,getNumberOfPartitions(tr));
  pInfo *partition = getPartition(chain, chain->frequRemem.model); 

  recordFrequRates(chain, chain->frequRemem.model, chain->frequRemem.numFrequRates, chain->frequRemem.curFrequRates);

  
  int state;
  double sum,curv;
  double r[chain->frequRemem.numFrequRates];
  
  chain->hastings=1;
  for(state = 0;state<chain->frequRemem.numFrequRates ; state ++)
    {
      curv = partition->frequencies[state];
      //r[state] =  drawRandDouble(); 
      r[state] =  drawRandBiUnif(chain,curv); 
    chain->hastings*=curv/r[state];
    }
    
  sum=0;
  
  for(state = 0;state<chain->frequRemem.numFrequRates ; state ++)
    {
      sum+=r[state]; 
    }
    for(state = 0;state<chain->frequRemem.numFrequRates ; state ++)
    {
      partition->frequencies[state]=r[state]/sum; 
    }
  //recalculate eigens

  exa_initReversibleGTR(chain, chain->frequRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

/* chain->curprior=get_frequency_prior(chain, tr->partitionData[chain->frequRemem.model].frequencies); */
/* chain->newprior=get_frequency_prior(chain, chain->frequRemem.curFrequRates); */

  evaluateOnePartition(chain, tr->start, TRUE, chain->frequRemem.model);
}


void frequency_proposal_reset(state * chain)
{
  restore_frequ_rates(chain, chain->frequRemem.adef, chain->frequRemem.model, chain->frequRemem.numFrequRates, chain->frequRemem.curFrequRates);
}



/*
 * should be sliding window proposal
 */

void edit_subs_rates(state *chain, int model, int subRatePos, double subRateValue)
{
  pInfo *partition = getPartition(chain,model); 

  assert(partition->dataType = DNA_DATA);
  assert(subRateValue <= RATE_MAX && subRateValue >= RATE_MIN);
  int states = partition->states; 
  int numSubsRates = (states * states - states) / 2;
  assert(subRatePos >= 0 && subRatePos < numSubsRates);
  partition->substRates[subRatePos] = subRateValue;
}



static void simple_model_proposal_reset(state * chain)
{
  restore_subs_rates(chain, chain->modelRemem.adef, chain->modelRemem.model, chain->modelRemem.numSubsRates, chain->modelRemem.curSubsRates);
}

nodeptr select_random_subtree(state *chain, tree *tr)
{
  nodeptr 
    p;

  do
    {
      int 
        exitDirection = drawRandInt(chain,3); 
     
      int r = drawRandInt(chain,tr->mxtips - 2) ; 

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



void printProposalType(proposal_type which_proposal)
{
  switch(which_proposal)
    {
    case E_SPR: 
      printf("E_SPR\n");
      break; 
    case E_SPR_MAPPED: 
      printf("E_SPR_MAPPED\n");
      break; 
    case UPDATE_MODEL : 
      printf("UPDATE_MODEL \n");
      break; 
    case UPDATE_GAMMA : 
      printf("UPDATE_GAMMA \n");
      break; 
    case UPDATE_GAMMA_EXP: 
      printf("UPDATE_GAMMA_EXP\n");
      break; 
    case UPDATE_SINGLE_BL: 
      printf("UPDATE_SINGLE_BL\n");
      break; 
    case UPDATE_SINGLE_BL_EXP : 
      printf("UPDATE_SINGLE_BL_EXP \n");
      break; 
    case UPDATE_SINGLE_BL_BIUNIF: 
      printf("UPDATE_SINGLE_BL_BIUNIF\n");
      break; 
    case UPDATE_MODEL_BIUNIF: 
      printf("UPDATE_MODEL_BIUNIF\n");
      break; 
    case UPDATE_MODEL_SINGLE_BIUNIF: 
      printf("UPDATE_MODEL_SINGLE_BIUNIF\n");
      break; 
    case UPDATE_MODEL_ALL_BIUNIF: 
      printf("UPDATE_MODEL_ALL_BIUNIF\n");
      break; 
    case UPDATE_FREQUENCIES_BIUNIF: 
      printf("UPDATE_FREQUENCIES_BIUNIF\n");
      break; 
    case UPDATE_MODEL_PERM_BIUNIF: 
      printf("UPDATE_MODEL_PERM_BIUNIF\n");
      break;  
//PROPOSALADD printProposalType NOTE Do not remove/modify  this line except for numerical value. The script addProposal.pl needs it as an identifier.
      

    default : 
      assert(0); 
    }
}





void getProposalFunctions(proposal_type ptype, proposal_functions* pF)
{
  /* TODO proposal add for script  */

 
  switch(ptype)
    {
    case  UPDATE_MODEL:
      pF->ptype =  UPDATE_MODEL;
      pF->pSubType = STANDARD; 
      pF->apply_func =  simple_model_proposal_apply;
      pF->reset_func =  simple_model_proposal_reset;
      pF->get_prior_ratio = get_branch_length_prior;
      break;

    case  UPDATE_SINGLE_BL:
      pF->ptype = UPDATE_SINGLE_BL;
      pF->pSubType = STANDARD; 
      pF->apply_func =  random_branch_length_proposal_apply;
      pF->reset_func =  random_branch_length_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior ;
      break;

    case  UPDATE_SINGLE_BL_EXP:
      pF->ptype = UPDATE_SINGLE_BL_EXP;
      pF->pSubType = EXP_DISTR; 
      pF->apply_func	=  random_branch_length_proposal_apply;
      pF->reset_func =  random_branch_length_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior ;
      break;

    case  UPDATE_GAMMA:
      pF->ptype = UPDATE_GAMMA;
      pF->pSubType = STANDARD; 
      pF->apply_func	=  simple_gamma_proposal_apply;
      pF->reset_func =  simple_gamma_proposal_reset;
      pF->get_prior_ratio =   get_alpha_prior;
      break;

    case  UPDATE_GAMMA_EXP:
      pF->ptype = UPDATE_GAMMA_EXP;
      pF->pSubType = EXP_DISTR; 
      pF->apply_func	=  simple_gamma_proposal_apply;
      pF->reset_func =  simple_gamma_proposal_reset;
      pF->get_prior_ratio =   get_alpha_prior;
      break;

    case  UPDATE_SINGLE_BL_BIUNIF:
      pF->ptype = UPDATE_SINGLE_BL_BIUNIF;
      pF->pSubType = BIUNIF_DISTR; 
      pF->apply_func	=  random_branch_length_proposal_apply;
      pF->reset_func =  random_branch_length_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior;
      break;

    case  UPDATE_MODEL_BIUNIF:
      pF->ptype = UPDATE_MODEL_BIUNIF;
      pF->pSubType = BIUNIF_DISTR; 
      pF->apply_func	=  simple_model_proposal_apply;
      pF->reset_func =  simple_model_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior;
      break;
	
    case  UPDATE_MODEL_SINGLE_BIUNIF:
      pF->ptype = UPDATE_MODEL_SINGLE_BIUNIF;
      pF->pSubType = STANDARD; 
      pF->apply_func	=  single_biunif_model_proposal_apply;
      pF->reset_func =  simple_model_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior;
      break;

    case  UPDATE_MODEL_ALL_BIUNIF:
      pF->ptype = UPDATE_MODEL_ALL_BIUNIF;
      pF->pSubType = STANDARD; 
      pF->apply_func	=  all_biunif_model_proposal_apply;
      pF->reset_func =  simple_model_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior;
      break;

    case  UPDATE_MODEL_PERM_BIUNIF:
      pF->ptype = UPDATE_MODEL_PERM_BIUNIF;
      //pF->pSubType = BIUNIF_PERM_DISTR; 
      pF->pSubType = STANDARD; 
      pF->apply_func	=  perm_biunif_model_proposal_apply;
      //pF->apply_func	=  simple_model_proposal_apply;
      pF->reset_func =  simple_model_proposal_reset;
      pF->get_prior_ratio =  get_branch_length_prior;
      break;
	
    case  UPDATE_FREQUENCIES_BIUNIF:
      pF->ptype = UPDATE_FREQUENCIES_BIUNIF;
      pF->pSubType = STANDARD; 
      pF->apply_func	=  frequency_proposal_apply;
      pF->reset_func =  frequency_proposal_reset;
      pF->get_prior_ratio =  get_frequency_prior;
      break;


    case  E_SPR:
      pF->ptype = E_SPR;      
      //pF->pSubType = STANDARD; //TODO change back to this and implement other option as propposal
      pF->pSubType = SPR_ADJUST;
      pF->apply_func =  extended_spr_apply;
      pF->reset_func =  extended_spr_reset;
      pF->get_prior_ratio =   get_branch_length_prior;
      break;

    case E_SPR_MAPPED:
      pF->ptype = E_SPR_MAPPED; 
      pF->pSubType = SPR_MAPPED; 
      pF->apply_func =  extended_spr_apply;
      pF->reset_func =  extended_spr_reset;
      pF->get_prior_ratio =   get_branch_length_prior;
      break;
      
      //PROPOSALADD getProposalFunctions NOTE Do not remove/modify  this line except for numerical value. The script addProposal.pl needs it as an identifier.

    default : 
      assert(0); 
    }
}


/* so here the idea would be to randomly choose among proposals? we can use typedef enum to label each, and return that */ 
static proposal_type select_proposal_type(state * chain)
{
  chain->newprior = chain->brLenRemem.bl_prior; //TODO Why is this here? 
  return drawSampleProportionally(chain,chain->proposalWeights, NUM_PROPOSALS) ; 
}




void resetSuccessCounter(int runid)
{
  state *start =  gAInfo.allChains + (runid * gAInfo.numberCoupledChains) ; 
  
  for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
    {
      state *chain = start + i ; 
      memset(chain->rejectedProposals, 0 , sizeof(int) * NUM_PROPOSALS); 
      memset(chain->acceptedProposals, 0 , sizeof(int) * NUM_PROPOSALS); 
    }
}


void step(state *chain)
{
  tree *tr = chain->tr;   

  proposal_type which_proposal;
  double testr;
  double acceptance;

  // just for validation (make sure we compare the same)
  evaluateGenericWrapper(chain, tr->start, FALSE);

  tr->startLH = tr->likelihood;

  // select proposal type
  which_proposal = select_proposal_type( chain );
    
  proposal_functions pF; 
  getProposalFunctions(which_proposal, &pF); 


  // apply the proposal function  
  pF.apply_func(chain, pF.pSubType);

  // FIXME: why is this here?
  if (chain->currentGeneration == 0 )
    {
      chain->curprior = chain->newprior;
    }
  //     PRINT("proposal done, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);

  //proposalTime += gettime() - t;
  /* decide upon acceptance */
  testr = drawRandDouble01(chain);

  //should look something like 
  /* acceptance = fmin(1,(chain->hastings) * 
     (exp(chain->newprior-chain->curprior)) * (exp(chain->tr->likelihood-chain->tr->startLH)));*/ 

  double myHeat = getChainHeat(chain ) ; 

  acceptance = fmin(1,(chain->hastings) * 
		    /* TODO for chain swapping ratio must be replaced again by proper prior   */
		    pF.get_prior_ratio(chain) 
		    /* (chain->newprior/chain->curprior) */
		    * (exp((tr->likelihood - tr->startLH) * myHeat)  ) 
		    );


  assert(which_proposal < NUM_PROPOSALS); 
      
  if(testr < acceptance)
    {
#ifdef DEBUG_SHOW_EACH_PROPOSAL
      if(processID == 0)
	{
	  printInfo(chain, "accepting\t");   
	  printProposalType(which_proposal); 
	}
#endif
      chain->acceptedProposals[which_proposal]++; 
      chain->totalAccepted++;
      // chain->proposalWeights[which_proposal] /= chain->penaltyFactor;
      penalize(chain, which_proposal, 1);

      tr->startLH = tr->likelihood;  //new LH
      chain->curprior = chain->newprior;          
    }
  else
    {
#ifdef DEBUG_SHOW_EACH_PROPOSAL
      if(processID == 0)
	{
	  printInfo(chain, "rejecting\t");   
	  printProposalType(which_proposal); 
	}
#endif
      pF.reset_func(chain); 
      /* dispatch_proposal_reset(which_proposal,chain); */
      chain->rejectedProposals[which_proposal]++;
      chain->totalRejected++;
      //chain->proposalWeights[which_proposal] *= chain->penaltyFactor; 
      penalize(chain, which_proposal, 0);



#if 0       
      /* probably a bad idea to comment that out. But let's assume
	 that we can do a verification, whenever we call evaluate. 
	 Then this here is not necessary.  */

      expensiveVerify(tr);

      // just for validation
      if(fabs(tr->startLH - tr->likelihood) > 1.0E-15)//TODO change back to 1.0E-15
      	{
      	  PRINT("WARNING: LH diff %.20f\n", tr->startLH - tr->likelihood);
      	  PRINT("after reset, iter %d tr LH %f, startLH %f\n", chain->currentGeneration, tr->likelihood, tr->startLH);
      	}
      assert(fabs(tr->startLH - tr->likelihood) < 0.1);
#endif
    }
  

#ifdef DEBUG_LNL_VERIFY
  int count = 0; 
  traverseAndCount(tr->start->back, &count, tr); 
  if(count != 2 * tr->mxtips - 3 )
    {      
      char tmp[10000]; 
      Tree2stringNexus(tmp, chain, tr->start->back, 0); 
      if(processID==0)
	printf("faulty TOPOLOGY: %s\n", tmp);

      assert(2 * tr->mxtips-3 == count); 
    }
#endif


  if(chain->couplingId == 0 && (chain->currentGeneration % chain->samplingFrequency) == chain->samplingFrequency - 1  )
    {
      if(processID == 0)
	{
	  printSample(chain);       

 	  chainInfo(chain); 
	  /* chainInfoOutput(chain);  // , sum_radius_accept, sum_radius_reject      	   */
	  resetSuccessCounter(chain->id / gAInfo.numberCoupledChains);
	}

      if(chain->currentGeneration >  BURNIN)
	addBipartitionsToHash(tr, chain ); 

    }

  chain->likelihood = tr->likelihood; 
  
  chain->currentGeneration++; 
}



/****************************************************************************************/
/* garbage : this is code that was not used in the initial version, but may be of use   */
/****************************************************************************************/





#if 0 


/*
 * should be sliding window proposal
 */

static boolean simpleBranchLengthProposalApply(state * chain)
{
   
  //for each branch get the current branch length
  //pull a uniform like
  //x = current, w =window
  //uniform(x-w/2,x+w/2)

  update_all_branches(chain, FALSE);
  evaluateGenericWrapper(chain, chain->tr->start, TRUE); /* update the tr->likelihood */

  //for prior, just using exponential for now
  //calculate for each branch length
  // where lambda is chosen and x is the branch length
  //lambda * exp(-lamba * x)

  //only calculate the new ones
  //
  return TRUE;
}

static void simpleBranchLengthProposalReset(state * chain)
{
  update_all_branches(chain, TRUE);
}

static void update_all_branches(state * s, boolean resetBL)
{
  int updated_branches = 0;
  assert(isTip(s->tr->start->number, s->tr->mxtips));
  /* visit each branch exactly once */
  traverse_branches(s->tr->start->back, &updated_branches, s, resetBL);
  assert(updated_branches == s->tr->mxtips + s->tr->mxtips - 3);
}


static void traverse_branches(nodeptr p, int *count, state * s, boolean resetBL)
{
  nodeptr q;
  //printf("current BL at %db%d: %f\n", p->number, p->back->number, p->z[0]);
  if(resetBL)
    reset_branch_length(p, s->tr->numBranches);
  else//can allow for other methods later
    set_branch_length_sliding_window(p, s->tr->numBranches, s, TRUE);
  *count += 1;


  if (! isTip(p->number, s->tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverse_branches(q->back, count, s, resetBL);
	  q = q->next;
	}   
      // WTF? in each recursion?
      newviewGeneric(s->tr, p, FALSE);     // not sure if we need this
    }
}


#endif
