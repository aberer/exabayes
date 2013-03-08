#include "common.h"

#include "axml.h"
#include "proposalStructs.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"

#include "convergence.h"


#include "chain.h"

#include "utils-topo.h" 


void expensiveVerify(tree *tr); 


/* static int spr_depth = 0; */

nodeptr select_random_subtree(tree *tr); 
void edit_subs_rates(tree *tr, int model, int subRatePos, double subRateValue);






void traverseAndCount(nodeptr p, int *count, tree *tr )
{
  nodeptr q;  
  assert(tr->numBranches > 0); 

  *count += 1;

  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseAndCount(q->back, count, tr);
	  q = q->next;
	} 
    }
}





//reads proposalWeights p and sets t for logistic function such that f(t)=p
void findLogisticT(state *curstate){
  static double leftShift=2.0;
   for(int i = 0; i < NUM_PROPOSALS;++i){
     if(curstate->proposalWeights[i]>0){
    curstate->proposalLogisticT[i]=-log((1.0-curstate->proposalWeights[i])/curstate->proposalWeights[i])+leftShift;  
    //printf("p(i): %f e^..: %f\n",curstate->proposalWeights[i], (1.0/(1.0+exp(leftShift-curstate->proposalLogisticT[i])) ));
     }else{
       curstate->proposalLogisticT[i]=-1.0;
     }
   }
}

//get values for logistic function
void findLogisticP(state *curstate){
  static double leftShift=2.0;
  for(int i = 0; i < NUM_PROPOSALS;++i){
    if(curstate->proposalLogisticT[i]>=0){
    curstate->proposalWeights[i]=(1.0/(1.0+exp(leftShift-curstate->proposalLogisticT[i])) );  
  //  printf("p(i): %f e^..: %f t: %f\n",curstate->proposalWeights[i], (1.0/(1.0+exp(leftShift-curstate->proposalLogisticT[i])) ), curstate->proposalLogisticT[i]);
   }
  }
}



static void recordSubsRates(tree *tr, int model, int numSubsRates, double *prevSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = tr->partitionData[model].substRates[i];
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

int extended_spr_traverse(state *curstate, nodeptr *insertNode)
{
  int 
    result, 
    r = drawRandDouble(); 
  
  radius++;   

  /* *insertNode = NULL;  */
  if(r < 0.5 )
    {


       *insertNode = isTip((*insertNode)->next->number, curstate->tr->mxtips) ? (*insertNode)->next :  (*insertNode)->next->back; 

     // *insertNode = (*insertNode)->next->back ;
    }
  else 
    {
       *insertNode = isTip((*insertNode)->next->next->number, curstate->tr->mxtips) ? (*insertNode)->next->next :  (*insertNode)->next->next->back; 
      //*insertNode = (*insertNode)->next->next->back;
    }

  /*
  if(isTip( (*insertNode)->number, curstate->tr->mxtips))
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
    
     if(isTip( (*insertNode)->number, curstate->tr->mxtips))
  {
     if(randNum == 0)
    direction--;
    else
    direction++;
  }
  */
  double randprop = drawRandDouble();
  result = randprop < curstate->eSprStopProb; 
  
  return result; 
}








static void extended_spr_apply(state *instate, int pSubType)
{
  tree * tr = instate->tr;

#ifdef DEBUG_SHOW_TREE
  char tmp[10000]; 
  Tree2stringNexus(tmp, tr, tr->start->back, 0); 
  if(processID==0)
    printf("topo before: %s\n", tmp);
#endif

  
  nodeptr    
    p = select_random_subtree(tr);

#if 0
  parsimonySPR(p, tr);
#endif  

  instate->sprMoveRemem.p = p;
  instate->sprMoveRemem.nb  = p->next->back;
  instate->sprMoveRemem.nnb = p->next->next->back;
  
  record_branch_info(instate->sprMoveRemem.nb, instate->sprMoveRemem.nbz, instate->tr->numBranches);
  record_branch_info(instate->sprMoveRemem.nnb, instate->sprMoveRemem.nnbz, instate->tr->numBranches);
  

  /* initial remapping of BL of nodes adjacent to pruned node  */
  double zqr[NUM_BRANCHES];

instate->hastings = 1; 
    
  for(int i = 0; i < tr->numBranches; i++)
    {
      if(pSubType == STANDARD)
	{
	  zqr[i] = instate->sprMoveRemem.nb->z[i] * instate->sprMoveRemem.nnb->z[i];
	
	  instate->hastings *= log(zqr[i]);
	}
      else if(pSubType == SPR_MAPPED)
	{
	  zqr[i] = instate->sprMoveRemem.nb->z[i] ; 
	}
      else if(pSubType == SPR_ADJUST)
	{
	  zqr[i] = instate->sprMoveRemem.nb->z[i] * instate->sprMoveRemem.nnb->z[i];  
	  instate->hastings *= log(zqr[i]);
	  zqr[i] = sqrt(zqr[i]);
	}
      if(zqr[i] > zmax) zqr[i] = zmax;
      if(zqr[i] < zmin) zqr[i] = zmin;
    }

  hookup(instate->sprMoveRemem.nb, instate->sprMoveRemem.nnb, zqr, tr->numBranches); 
  p->next->next->back = p->next->back = (node *) NULL;
  /* done remove node p (omitted BL opt) */

  nodeptr initNode = NULL; 
  nodeptr curNode = NULL; 
  boolean remapBL = FALSE; 
  if(drawRandDouble() > 0.5)
    {
      curNode = instate->sprMoveRemem.nb; 
      remapBL = TRUE ; 
    }
  else 
    {
      curNode = instate->sprMoveRemem.nnb; 
    }
  initNode = curNode; 

  int accepted = FALSE;   

  radius = 0; 
//printf("curNode:  back next nextnext\n");
  while( NOT  accepted)
    {
      accepted = extended_spr_traverse(instate, &curNode); 
  //    if(processID == 0) 
//	printf("%d:  %d  %d  %d\n",curNode->number, curNode->back->number, curNode->next->back->number, curNode->next->next->back->number); 

      /* needed for spr remap */
      if(curNode == initNode)
	remapBL = NOT remapBL; 
    }



  instate->newprior = 1; 
  instate->curprior = 1; 

  instate->sprMoveRemem.q = curNode;
  instate->sprMoveRemem.r = instate->sprMoveRemem.q->back;
  record_branch_info(instate->sprMoveRemem.q, instate->brLenRemem.qz, instate->tr->numBranches);

  
  
   if(pSubType == STANDARD){
  for(int branchCount=0; branchCount<instate->tr->numBranches; branchCount++)
     {
      instate->hastings/=log(curNode->z[branchCount]);
     }
   }else if(pSubType == SPR_ADJUST){
  for(int branchCount=0; branchCount<instate->tr->numBranches; branchCount++)
     {
      instate->hastings/=(2*log(curNode->z[branchCount]));
     }
   }
  
 // Thorough = 0;
  // assert(tr->thoroughInsertion == 0);
  /* insertBIG wont change the BL if we are not in thorough mode */

  if(pSubType == STANDARD)
    {
     insertWithUnifBL(instate->sprMoveRemem.p, instate->sprMoveRemem.q, instate->tr->numBranches);
  /*   for(int branchCount=0; branchCount<instate->tr->numBranches; branchCount++)
     {
       instate->hastings/=log(instate->brLenRemem.qz[branchCount]);
     }
     */
   //  printf("hastings: %f\n", instate->hastings);
      // insertBIG(instate->tr, instate->sprMoveRemem.p, instate->sprMoveRemem.q, instate->tr->numBranches);
    }
  else if(pSubType == SPR_MAPPED)
    {
      double *neighborZ = remapBL ? instate->sprMoveRemem.nbz :  instate->sprMoveRemem.nnbz; 
      
        if( remapBL ) 
    {
      for(int i = 0; i < tr->numBranches; ++i)
	instate->sprMoveRemem.nb->z[i] = instate->sprMoveRemem.nb->back->z[i] = instate->sprMoveRemem.nnbz[i]; 
    }
      
      insertWithGenericBL(instate->sprMoveRemem.p, instate->sprMoveRemem.q, instate->sprMoveRemem.p->z, curNode->z, neighborZ, tr->numBranches);


      /* IMPORTANT TODO verify, that the mapping actually works, as we
	 had this in mind. Use topo print functions for that.
       */


      /* TODO PERFORMANCE */
      /* newviewGeneric(tr, instate->sprMoveRemem.p, FALSE); */
    }
     else if(pSubType == SPR_ADJUST)
    {
      insertWithUnifBLScaled(instate->sprMoveRemem.p, instate->sprMoveRemem.q, 2.0,  instate->tr->numBranches);
    }

//   if( (pSubType == SPR_MAPPED ) && remapBL) 
//     {
//       for(int i = 0; i < tr->numBranches; ++i)
// 	instate->sprMoveRemem.nb->z[i] = instate->sprMoveRemem.nb->back->z[i] = instate->sprMoveRemem.nnbz[i]; 
//     }

  /* TODO problem here? is not tip?? FIXED: should never be a tip anymore. Possible since we are interested in edges and each edge connects to at least one inner node*/
#if 0 
  evaluateGeneric(instate->tr, instate->sprMoveRemem.p->next->next, FALSE);
#else 
  evaluateGeneric(instate->tr, instate->tr->start, TRUE);
#endif
  expensiveVerify(tr); 

}





static void extended_spr_reset(state * instate)
{
  tree *tr = instate->tr; 

  /* prune the insertion */
  hookup(instate->sprMoveRemem.q, instate->sprMoveRemem.r, instate->brLenRemem.qz, instate->tr->numBranches);

  instate->sprMoveRemem.p->next->next->back = instate->sprMoveRemem.p->next->back = (nodeptr) NULL;

  /* insert the pruned tree in its original node */
  hookup(instate->sprMoveRemem.p->next,        instate->sprMoveRemem.nb, instate->sprMoveRemem.nbz, instate->tr->numBranches);
  hookup(instate->sprMoveRemem.p->next->next, instate->sprMoveRemem.nnb, instate->sprMoveRemem.nnbz, instate->tr->numBranches);
  
  if(processID == 0)
    {

#ifdef DEBUG_SHOW_TREE
      char tmp[100000];
      Tree2stringNexus(tmp, instate->tr, instate->tr->start->back, 0); 
      printf("topo reset: %s\n", tmp); 
#endif
    }

  evaluateGeneric(tr, tr->start, TRUE);

  newviewGeneric(instate->tr, instate->sprMoveRemem.p, FALSE);
  double val1 = instate->tr->likelihood; 
  
  newviewGeneric(instate->tr, instate->tr->start, TRUE);
  double  val2 = instate->tr->likelihood; 

  assert( fabs ( val2 - val1 ) < 0.0001 ); 

}


//--------Alpha-Proposal-for-GAMMA-----------------------------------------------------------------

double get_alpha_prior(state *curstate )
{
  
 return 1;//TODO obviously needs acctual prior 
}





static void simple_gamma_proposal_apply(state * instate, int pSubType)
{
  //TODO: add safety to max and min values
  double newalpha, curv, r,mx,mn;
  instate->modelRemem.model=drawRandInt(instate->tr->NumberOfModels);
  curv = instate->tr->partitionData[instate->modelRemem.model].alpha;
  instate->gammaRemem.curAlpha = curv;
  
  switch(pSubType)
        {
        case STANDARD://simple sliding window
	  r = drawRandDouble();
	  mn = curv-(instate->gammaRemem.gm_sliding_window_w/2);
	  mx = curv+(instate->gammaRemem.gm_sliding_window_w/2);
	  newalpha = fabs(mn + r * (mx-mn));
	  break;
        case EXP_DISTR:
          newalpha  = drawRandExp(1/curv);
          break;
        default:
          assert(0);
        }
  /* Ensure always you stay within this range */
  if(newalpha > ALPHA_MAX) newalpha = ALPHA_MAX;
  if(newalpha < ALPHA_MIN) newalpha = ALPHA_MIN;
  
  instate->hastings = 1; //since it is symmetrical, hastings=1
  instate->newprior = get_alpha_prior(instate); 
  instate->curprior = get_alpha_prior(instate); 
  
  instate->tr->partitionData[instate->modelRemem.model].alpha = newalpha;

  makeGammaCats(instate->tr->partitionData[instate->modelRemem.model].alpha, instate->tr->partitionData[instate->modelRemem.model].gammaRates, 4, instate->tr->useMedian);

  
  evaluateGeneric(instate->tr, instate->tr->start, TRUE);
}


// static void exp_gamma_proposal_apply(state * instate, int pSubType)
// {
//   double newalpha, curv;
//   curv = instate->tr->partitionData[instate->modelRemem.model].alpha;
//   instate->gammaRemem.curAlpha = curv;
//   newalpha  = drawRandExp(1/curv);
//   
// 
//  
//   /* Ensure always you stay within this range */
//   while(newalpha > ALPHA_MAX || newalpha < ALPHA_MIN){
//   if(newalpha > ALPHA_MAX) newalpha = 2*ALPHA_MAX-newalpha;
//   if(newalpha < ALPHA_MIN) newalpha = 2*ALPHA_MIN-newalpha;
//   }
//   instate->hastings = (1/newalpha)*exp(-(1/newalpha)*curv)/((1/curv)*exp(-(1/curv)*newalpha)); //TODO do not ignore reflection
//   /* instate->newprior = get_alpha_prior(instate);  */
//   /* instate->curprior = get_alpha_prior(instate);  */
//   
//   instate->tr->partitionData[instate->modelRemem.model].alpha = newalpha;
//   
//   makeGammaCats(instate->tr->partitionData[instate->modelRemem.model].alpha, instate->tr->partitionData[instate->modelRemem.model].gammaRates, 4, instate->tr->useMedian);
// 
//   evaluateGeneric(instate->tr, instate->tr->start, TRUE);
// }


static void simple_gamma_proposal_reset(state * instate)
{
  instate->tr->partitionData[instate->modelRemem.model].alpha = instate->gammaRemem.curAlpha;

  makeGammaCats(instate->tr->partitionData[instate->modelRemem.model].alpha, instate->tr->partitionData[instate->modelRemem.model].gammaRates, 4, instate->tr->useMedian);

  evaluateGeneric(instate->tr, instate->tr->start, TRUE);
}

//------------------------------------------------------------------------------

/* TODO we do not even use this function, do we? NOTE Now we do ;) */
void penalize(state *curstate, int which_proposal, int acceptance)
{
  double max=4.0;
  double min=0.0;
  if(acceptance){
    curstate->proposalLogisticT[which_proposal]+=curstate->penaltyFactor;
    if(curstate->proposalLogisticT[which_proposal]>max)
      curstate->proposalLogisticT[which_proposal]=max;
    //curstate->proposalWeights[which_proposal] /= curstate->penaltyFactor; 
  }
  else {
  /*  if(curstate->totalRejected>0){
    curstate->proposalLogisticT[which_proposal]-=((curstate->totalAccepted/curstate->totalRejected)*curstate->penaltyFactor);*///TODO check whether "VCG" is better than considering all
    if(curstate->totalRejected - curstate->rejectedProposals[which_proposal]>0){
    curstate->proposalLogisticT[which_proposal]-=(((curstate->totalAccepted- curstate->acceptedProposals[which_proposal])/(curstate->totalRejected - curstate->rejectedProposals[which_proposal]))*curstate->penaltyFactor);      
    }else{
      curstate->proposalLogisticT[which_proposal]-=curstate->penaltyFactor;
    }
    if(curstate->proposalLogisticT[which_proposal]<min)
      curstate->proposalLogisticT[which_proposal]=min;
    //curstate->proposalWeights[which_proposal] *= curstate->penaltyFactor; 
  }
  findLogisticP(curstate);
  //normalizeProposalWeights(curstate); 
}





//NOTE: should only be called at the very beginning. Afterwards a probability sum of 1.0 is not required.

void normalizeProposalWeights(state *curstate)
{
  double sum = 0 ; 
  for(int i = 0; i < NUM_PROPOSALS;++i)
    sum += curstate->proposalWeights[i]; 
  
  for(int i = 0; i < NUM_PROPOSALS;++i)
    curstate->proposalWeights[i] /= sum; 
  
  findLogisticT(curstate);
}



static void simple_model_proposal_apply(state *instate, int pSubType)//llpqr
{
  //TODO: add safety to max and min values
  //record the old ones
  instate->modelRemem.model=drawRandInt(instate->tr->NumberOfModels);
  recordSubsRates(instate->tr, instate->modelRemem.model, instate->modelRemem.numSubsRates, instate->modelRemem.curSubsRates);
  //choose a random set of model params,
  //probably with dirichlet proposal
  //with uniform probabilities, no need to have other
  int state, changeState;
  double new_value,curv;
  double r;
  
  

  
  double mx,mn; //for sliding window
  
  int list[instate->modelRemem.numSubsRates];//for biunif_distr and biunif_perm_distr
  
  int numberOfEdits;//for biunif_perm_distr
  
  switch(pSubType)
        {
	  case STANDARD:
	    case BIUNIF_DISTR:
	      numberOfEdits=instate->modelRemem.numSubsRates;
	      break;
	      case BIUNIF_PERM_DISTR:
	    numberOfEdits=drawRandInt(instate->modelRemem.numSubsRates);
		break;
	} 
  
  
  instate->hastings=1.0;
  
    for(state = 0;state<instate->modelRemem.numSubsRates ; state ++)
    {
      
      switch(pSubType)
        {
        case STANDARD: //using the branch length sliding window for a test    
	  changeState=state;
	  curv = instate->tr->partitionData[instate->modelRemem.model].substRates[state];
	  r =  drawRandDouble();
	  mn = curv-(instate->modelRemem.rt_sliding_window_w/2);
	  mx = curv+(instate->modelRemem.rt_sliding_window_w/2);
	
	  new_value = fabs(mn + r * (mx-mn));
	  
	  /* Ensure always you stay within this range */
      if(new_value > RATE_MAX) new_value = RATE_MAX;
      if(new_value < RATE_MIN) new_value = RATE_MIN;
      
	  break;
        
	case BIUNIF_DISTR:
	changeState=drawRandInt(instate->modelRemem.numSubsRates);
      if(list[changeState]!=1)
      {
      list[changeState]=1;;      
      curv = instate->tr->partitionData[instate->modelRemem.model].substRates[changeState];
      r =  drawRandBiUnif(curv);
      new_value = r;
      instate->hastings*=curv/new_value;
 
      while(new_value> RATE_MAX|| new_value< RATE_MIN)
      {
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      }
      break;
      case BIUNIF_PERM_DISTR:
	  drawPermutation(list, instate->modelRemem.numSubsRates);
	 curv = instate->tr->partitionData[instate->modelRemem.model].substRates[list[state]];
	 r =  drawRandBiUnif(curv);
       new_value = r;
      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       instate->hastings*=curv/new_value;
	 break;
	 
	case SINGLE_BIUNIF:
	  drawPermutation(list, instate->modelRemem.numSubsRates);
	 curv = instate->tr->partitionData[instate->modelRemem.model].substRates[list[state]];
	 r =  drawRandBiUnif(curv);
       new_value = r;
      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       instate->hastings*=curv/new_value;
	 break;
	 
        default:
          assert(0);
        
      

    }//end of switch

edit_subs_rates(instate->tr,instate->modelRemem.model, changeState, new_value);
    }
  //recalculate eigens

  initReversibleGTR(instate->tr, instate->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* TODO: need to broadcast rates here for parallel version ! */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
  //TODO: without this, the run will fail after a successful model, but failing SPR
  //TODOFER: what did we have in mind regarding the comment above?
  
  /* evaluateGeneric(instate->tr, instate->tr->start, FALSE); */
  //for prior, just use dirichlet
  // independent gamma distribution for each parameter
  //the pdf for this is
  // for gamma the prior is gamma

  //for statefreqs should all be uniform

  //only calculate the new ones
}


// //draws a random subset (drawing with replacement) of the states and changes the according to biunif distribution.
// static void biunif_model_proposal_apply(state *instate, int pSubType)
// {
//   //record the old one 
//    instate->modelRemem.model=drawRandInt(instate->tr->NumberOfModels);
//   recordSubsRates(instate->tr, instate->modelRemem.model, instate->modelRemem.numSubsRates, instate->modelRemem.curSubsRates);
//   int state, randState;
//   double new_value,curv;
//   double r;
//   int list[instate->modelRemem.numSubsRates];
//   
//   
//    instate->hastings=1.0;
//   for(state = 0;state<instate->modelRemem.numSubsRates ; state ++)
//     {
//       randState=drawRandInt(instate->modelRemem.numSubsRates);
//       if(list[randState]!=1)
//       {
//       list[randState]=1;;
//       
//       curv = instate->tr->partitionData[instate->modelRemem.model].substRates[randState];
//       r =  drawRandBiUnif(curv);
// 
//       new_value = r;
// 
//       while(new_value> RATE_MAX|| new_value< RATE_MIN){
//       if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
//       if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
//       }
//       
//        instate->hastings*=curv/new_value;
//       edit_subs_rates(instate->tr,instate->modelRemem.model, randState, new_value);
//       }
//     }
//       
// 
// 
//   initReversibleGTR(instate->tr, instate->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
// 
//   evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
// 
// }


// static void perm_biunif_model_proposal_apply(state *instate, int pSubType)
// {
//   //record the old one 
//    instate->modelRemem.model=drawRandInt(instate->tr->NumberOfModels);
//   recordSubsRates(instate->tr, instate->modelRemem.model, instate->modelRemem.numSubsRates, instate->modelRemem.curSubsRates);
//   int state, randNumber;
//   double new_value,curv;
//   double r;
//   
//   randNumber=drawRandInt(instate->modelRemem.numSubsRates);
//     int perm[instate->modelRemem.numSubsRates];
//   drawPermutation(perm, instate->modelRemem.numSubsRates);
//   
//    instate->hastings=1.0;
//   for(state = 0;state<randNumber ; state ++)
//     {
//       
//       
//       curv = instate->tr->partitionData[instate->modelRemem.model].substRates[perm[state]];
//       r =  drawRandBiUnif(curv);
// 
//       new_value = r;
// 
//       while(new_value> RATE_MAX|| new_value< RATE_MIN){
//       if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
//       if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
//       }
//       
//        instate->hastings*=curv/new_value;
//       edit_subs_rates(instate->tr,instate->modelRemem.model, perm[state], new_value);
//       
//     }
//       
//   
// 
//   initReversibleGTR(instate->tr, instate->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
// 
//   evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
//   
// }


static void single_biunif_model_proposal_apply(state *instate,int pSubType)//NOTE whenever a model parameter changes, all branch lengths have to be re-normalized with 1/fracchange. Additionally we always must do a full tree traversal to get the likelihood. So updating a single parameter is rather expensive, .
{
  //record the old one //TODO sufficient to store single value.
   instate->modelRemem.model=drawRandInt(instate->tr->NumberOfModels);
  recordSubsRates(instate->tr, instate->modelRemem.model, instate->modelRemem.numSubsRates, instate->modelRemem.curSubsRates);
  //choose a random set parameter,
  //with uniform probabilities

  int  randState=drawRandInt(instate->modelRemem.numSubsRates);

  double new_value,curv;
  double r;
  
  //int state=drawRandInt(instate->modelRemem.numSubsRates);
  

  curv = instate->tr->partitionData[instate->modelRemem.model].substRates[randState];
  r =  drawRandBiUnif(curv);

  new_value = r;
      
  while(new_value> RATE_MAX|| new_value< RATE_MIN){
    if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
    if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
  }

  edit_subs_rates(instate->tr,instate->modelRemem.model, randState, new_value);

  instate->hastings=curv/new_value;

  initReversibleGTR(instate->tr, instate->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
}

static void all_biunif_model_proposal_apply(state *instate, int pSubType)
{
  //record the old one 
   instate->modelRemem.model=drawRandInt(instate->tr->NumberOfModels);
  recordSubsRates(instate->tr, instate->modelRemem.model, instate->modelRemem.numSubsRates, instate->modelRemem.curSubsRates);
  //choose a random set parameter,
  //with uniform probabilities
  int state;
  double new_value,curv;
  double r;
  
  
   instate->hastings=1.0;

  for(state = 0;state<instate->modelRemem.numSubsRates ; state ++)
    {
      curv = instate->tr->partitionData[instate->modelRemem.model].substRates[state];
      r =  drawRandBiUnif(curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       instate->hastings*=curv/new_value;
      edit_subs_rates(instate->tr,instate->modelRemem.model, state, new_value);
    }
      
  

  initReversibleGTR(instate->tr, instate->modelRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */

}

static void restore_subs_rates(tree *tr, analdef *adef, int model, int numSubsRates, double *prevSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    tr->partitionData[model].substRates[i] = prevSubsRates[i];

  initReversibleGTR(tr, model);

  /* TODO need to broadcast rates here for parallel version */

  evaluateGeneric(tr, tr->start, TRUE);
}


//--------Branch-Length_Proposals---------------------------------------------------

double get_branch_length_prior( state *curstate)
{//TODO decide on sensible prior
  return 1;  
}

//setting this out to allow for other types of setting
static void set_branch_length_sliding_window(nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
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
      r = drawRandDouble();

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

static void set_branch_length_biunif(nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
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
      r = drawRandBiUnif(real_z);//draw from [real_z/2,2*real_z]

    
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

static void set_branch_length_exp(nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
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
	r = drawRandExp(1.0/real_z);
  
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




static void random_branch_length_proposal_apply(state * instate, int pSubType)
{

  const int num_branches = (instate->tr->mxtips * 2) - 3;
  int target_branch = drawRandInt(num_branches); 
  node *p = select_branch_by_id_dfs( instate->tr->start, target_branch, instate );
  
    switch(pSubType)
        {
        case STANDARD://simple sliding window
	     //for each branch get the current branch length
	    //pull a uniform like
	    //x = current, w =window
	    //uniform(x-w/2,x+w/2)
	    set_branch_length_sliding_window(p, instate->tr->numBranches, instate, TRUE);
	  break;
        case EXP_DISTR:
	  set_branch_length_exp(p, instate->tr->numBranches, instate, TRUE);
          break;
	case BIUNIF_DISTR:
	set_branch_length_biunif(p, instate->tr->numBranches, instate, TRUE);
	  break;
        default:
          assert(0);
        }


  instate->brLenRemem.single_bl_branch = target_branch;
  evaluateGeneric(instate->tr, p, FALSE); 
  //evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* update the tr->likelihood *//FALSE seems to work
}


// static void biunif_branch_length_proposal_apply(state * instate, int pSubType)
// {
//    
//   //for one branch get the current branch length
//   //pull a uniform like
//   //x = current,
//   //uniform(x/2,x*2)
//   
//   
//   const int num_branches = (instate->tr->mxtips * 2) - 3;
//   int target_branch = drawRandInt(num_branches); 
//   node *p = select_branch_by_id_dfs( instate->tr->start, target_branch, instate );
//  
//   //set_branch_length_sliding_window(p, instate->tr->numBranches, instate, TRUE);
//   set_branch_length_biunif(p, instate->tr->numBranches, instate, TRUE);
// 
//   instate->brLenRemem.single_bl_branch = target_branch;
//   evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* update the tr->likelihood *///TODO see below
//   //   return TRUE;
// }

// static void exp_branch_length_proposal_apply(state * instate, int pSubType)
// {
//   const int num_branches = (instate->tr->mxtips * 2) - 3;
//   int target_branch = drawRandInt(num_branches); 
//   node *p = select_branch_by_id_dfs( instate->tr->start, target_branch, instate );
//   
//   set_branch_length_exp(p, instate->tr->numBranches, instate, TRUE);
// 
//   instate->brLenRemem.single_bl_branch = target_branch;
//   evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* update the tr->likelihood *///TODO see below
// }

static void random_branch_length_proposal_reset(state * instate)
{
  node *p;
  assert( instate->brLenRemem.single_bl_branch != -1 );
  
  // ok, maybe it would be smarter to store the node ptr for rollback rather than re-search it...
  p = select_branch_by_id_dfs( instate->tr->start, instate->brLenRemem.single_bl_branch, instate );
  
  reset_branch_length(p, instate->tr->numBranches);
  //   printf( "reset bl: %p %f\n", p, p->z[0] );
  //update_all_branches(instate, TRUE);
  evaluateGeneric(instate->tr, instate->tr->start, FALSE);
 // evaluateGeneric(instate->tr, p, FALSE); //This yields a very slight likelihood difference.NOTE if we want exact likelihoods as before the proposal, we must evaluate from instate->tr->start, that is: evaluateGeneric(instate->tr, instate->tr->start, TRUE);
  instate->brLenRemem.single_bl_branch = -1;
}

double get_frequency_prior(state * instate)
{
 return 1; 
}

static void restore_frequ_rates(tree *tr, analdef *adef, int model, int numFrequRates, double *prevFrequRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    tr->partitionData[model].frequencies[i] = prevFrequRates[i];

  initReversibleGTR(tr, model);

  evaluateGeneric(tr, tr->start, TRUE);
}

static void recordFrequRates(tree *tr, int model, int numFrequRates, double *prevFrequRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    prevFrequRates[i] = tr->partitionData[model].frequencies[i];
}

void frequency_proposal_apply(state * instate, int pSubType)
{
   instate->frequRemem.model=drawRandInt(instate->tr->NumberOfModels);
  recordFrequRates(instate->tr, instate->frequRemem.model, instate->frequRemem.numFrequRates, instate->frequRemem.curFrequRates);

  
  int state;
  double sum,curv;
  double r[instate->frequRemem.numFrequRates];
  
  instate->hastings=1;
  for(state = 0;state<instate->frequRemem.numFrequRates ; state ++)
    {
      curv = instate->tr->partitionData[instate->frequRemem.model].frequencies[state];
      //r[state] =  drawRandDouble(); 
    r[state] =  drawRandBiUnif(curv); 
    instate->hastings*=curv/r[state];
    }
    
  sum=0;
  
  for(state = 0;state<instate->frequRemem.numFrequRates ; state ++)
    {
      sum+=r[state]; 
    }
    for(state = 0;state<instate->frequRemem.numFrequRates ; state ++)
    {
      instate->tr->partitionData[instate->frequRemem.model].frequencies[state]=r[state]/sum; 
    }
  //recalculate eigens

  initReversibleGTR(instate->tr, instate->frequRemem.model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

/* instate->curprior=get_frequency_prior(instate, instate->tr->partitionData[instate->frequRemem.model].frequencies); */
/* instate->newprior=get_frequency_prior(instate, instate->frequRemem.curFrequRates); */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */

  //evaluateGeneric(instate->tr, instate->tr->start, FALSE);

}

void frequency_proposal_reset(state * instate)
{
  restore_frequ_rates(instate->tr, instate->frequRemem.adef, instate->frequRemem.model, instate->frequRemem.numFrequRates, instate->frequRemem.curFrequRates);
}



/*
 * should be sliding window proposal
 */

void edit_subs_rates(tree *tr, int model, int subRatePos, double subRateValue)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  assert(subRateValue <= RATE_MAX && subRateValue >= RATE_MIN);
  int states = tr->partitionData[model].states; 
  int numSubsRates = (states * states - states) / 2;
  assert(subRatePos >= 0 && subRatePos < numSubsRates);
  tr->partitionData[model].substRates[subRatePos] = subRateValue;
}



static void simple_model_proposal_reset(state * instate)
{
  restore_subs_rates(instate->tr, instate->modelRemem.adef, instate->modelRemem.model, instate->modelRemem.numSubsRates, instate->modelRemem.curSubsRates);
  //evaluateGeneric(instate->tr, instate->tr->start, FALSE);
}

nodeptr select_random_subtree(tree *tr)
{
  nodeptr 
    p;

  do
    {
      int 
        exitDirection = drawRandInt(3); 
     
      int r = drawRandInt(tr->mxtips - 2) ; 

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

      /* TODO proposal add PROPOSALADD  */

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
      pF->pSubType = SINGLE_BIUNIF; 
      pF->apply_func	=  simple_model_proposal_apply;
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
      pF->pSubType = BIUNIF_PERM_DISTR; 
      pF->apply_func	=  simple_model_proposal_apply;
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

    default : 
      assert(0); 
    }
}


/* so here the idea would be to randomly choose among proposals? we can use typedef enum to label each, and return that */ 
static proposal_type select_proposal_type(state * instate)
{
  instate->newprior = instate->brLenRemem.bl_prior; //TODO Why is this here? 
  return drawSampleProportionally(instate->proposalWeights, NUM_PROPOSALS) ; 
}



void step(state *curstate)
{
  tree *tr = curstate->tr;   

  proposal_type which_proposal;
  /* double t = gettime();  */
  /* double proposalTime = 0.0; */
  double testr;
  double acceptance;

  // just for validation (make sure we compare the same)
  evaluateGeneric(tr, tr->start, FALSE);
  expensiveVerify(tr); 

  tr->startLH = tr->likelihood;

  // select proposal type
  which_proposal = select_proposal_type( curstate );
    
  proposal_functions pF; 
  getProposalFunctions(which_proposal, &pF); 


  // apply the proposal function  
  pF.apply_func(curstate, pF.pSubType);

  // FIXME: why is this here?
  if (curstate->currentGeneration == 0 )
    {
      curstate->curprior = curstate->newprior;
    }
  //     PRINT("proposal done, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);

  //proposalTime += gettime() - t;
  /* decide upon acceptance */
  testr = drawRandDouble();

  //should look something like 
  /* acceptance = fmin(1,(curstate->hastings) * 
     (exp(curstate->newprior-curstate->curprior)) * (exp(curstate->tr->likelihood-curstate->tr->startLH)));*/

  acceptance = fmin(1,(curstate->hastings) * 
		    pF.get_prior_ratio(curstate) 
		    /* (curstate->newprior/curstate->curprior) */
		    * (exp(tr->likelihood - tr->startLH)));

  assert(which_proposal < NUM_PROPOSALS); 
      
  if(testr < acceptance)
    {
#ifdef DEBUG_SHOW_EACH_PROPOSAL
      if(processID == 0)
	{
	  printInfo(curstate, "accepting\t");   
	  printProposalType(which_proposal); 
	}
#endif
      curstate->acceptedProposals[which_proposal]++; 
      curstate->totalAccepted++;
      // curstate->proposalWeights[which_proposal] /= curstate->penaltyFactor;
      penalize(curstate, which_proposal, 1);

      tr->startLH = tr->likelihood;  //new LH
      curstate->curprior = curstate->newprior;          
    }
  else
    {
#ifdef DEBUG_SHOW_EACH_PROPOSAL
      if(processID == 0)
	{
	  printInfo(curstate, "rejecting\t");   
	  printProposalType(which_proposal); 
	}
#endif
      pF.reset_func(curstate); 
      /* dispatch_proposal_reset(which_proposal,curstate); */
      curstate->rejectedProposals[which_proposal]++;
      curstate->totalRejected++;
      //curstate->proposalWeights[which_proposal] *= curstate->penaltyFactor; 
      penalize(curstate, which_proposal, 0);
      
      /* TODO remove, if not needed any more */
/* #if 0  */
      /* evaluateGeneric(tr, tr->start, FALSE); */
/* #else  */
/* #endif */
  
      expensiveVerify(tr); 

      // just for validation 
      if(fabs(tr->startLH - tr->likelihood) > 1.0E-15)//TODO change back to 1.0E-15
	{
	  PRINT("WARNING: LH diff %.20f\n", tr->startLH - tr->likelihood);
	  PRINT("after reset, iter %d tr LH %f, startLH %f\n", curstate->currentGeneration, tr->likelihood, tr->startLH);
	}      
      assert(fabs(tr->startLH - tr->likelihood) < 0.1);
    }
  

#ifdef DEBUG_LNL_VERIFY
  int count = 0; 
  traverseAndCount(tr->start->back, &count, tr); 
  if(count != 2 * tr->mxtips - 3 )
    {      
      char tmp[10000]; 
      Tree2stringNexus(tmp, tr, tr->start->back, 0); 
      if(processID==0)
	printf("faulty TOPOLOGY: %s\n", tmp);

      assert(2 * tr->mxtips-3 == count); 
    }
#endif


  if((curstate->currentGeneration % curstate->samplingFrequency) == 0)
    {
      if(processID == 0)
	{
	  printSample(curstate);       
	  chainInfoOutput(curstate);  // , sum_radius_accept, sum_radius_reject      	  
	}
      addBipartitionsToHash(tr, curstate ); 
    }
  
  curstate->currentGeneration++; 
}



/****************************************************************************************/
/* garbage : this is code that was not used in the initial version, but may be of use   */
/****************************************************************************************/





#if 0 


/*
 * should be sliding window proposal
 */

static boolean simpleBranchLengthProposalApply(state * instate)
{
   
  //for each branch get the current branch length
  //pull a uniform like
  //x = current, w =window
  //uniform(x-w/2,x+w/2)

  update_all_branches(instate, FALSE);
  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* update the tr->likelihood */

  //for prior, just using exponential for now
  //calculate for each branch length
  // where lambda is chosen and x is the branch length
  //lambda * exp(-lamba * x)

  //only calculate the new ones
  //
  return TRUE;
}

static void simpleBranchLengthProposalReset(state * instate)
{
  update_all_branches(instate, TRUE);
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
