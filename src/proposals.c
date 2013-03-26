/**
   @file proposals.c

   @brief All proposals that make the chains of ExaBayes move in its
   space.
    
 */ 

#include "axml.h"
#include "bayes.h"
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
#include "misc-utils.h"
#include "branch.h"

void expensiveVerify(tree *tr); 

/* nodeptr select_random_subtree(state *chain, tree *tr); */
void edit_subs_rates(state *chain, int model, int subRatePos, double subRateValue);


#if 0 

/* quarantining this, until change is over */

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
#endif


static void recordSubsRates(state *chain, int model, int numSubsRates, double *prevSubsRates)
{
  pInfo *partition = getPartition(chain,model); 
  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = partition->substRates[i];
}


static void reset_branch_length(nodeptr p, int numBranches, double *bls)
{
  int i;
  for(i = 0; i < numBranches; i++)
    p->z[i] = p->back->z[i] = bls[i];   /* restore saved value */
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




/* TODO more diagnostics on where the stuff moves to  */
int extended_spr_traverse(state *chain, nodeptr *insertNode, double stopProp)
{
  int 
    result, 
    r = drawRandDouble01(chain); 

  if(r < 0.5 )
    *insertNode = (*insertNode)->next->back; 
  else 
    *insertNode = (*insertNode)->next->next->back; 

  double randprop = drawRandDouble01(chain);

  result = randprop < stopProp;
  
  return result; 
}



static nodeptr getThirdNode(tree *tr, int node, int neighBourA, int neighBourB)
{
  nodeptr 
    target = tr->nodep[node]; 

  if(target->back->number != neighBourA && target->back->number != neighBourB)
    return target; 
  else if (target->next->back->number != neighBourA && target->next->back->number != neighBourB)
    return target->next; 
  else 
    {
      assert(target->next->next->back->number != neighBourB && target->next->next->back->number != neighBourA); 
      return target->next->next; 
    }
  return NULL; 
}




static void onePartitionEval(state *chain, proposalFunction *thisProposal)
{
  int model = thisProposal->remembrance.partInfo->modelNum;
  evaluateOnePartition(chain, chain->tr->start, TRUE ,model);
}


static void dummy_eval(state *chain,proposalFunction *thisProposal)
{
  evaluateGenericWrapper(chain, chain->tr->start, TRUE);
}





static void spr_eval(state *chain, proposalFunction *thisProposal)
{  

#if 1 
  evaluateGenericWrapper(chain, chain->tr->start, TRUE); /* TODO   */
#else 
  /* TODO still so much to do */
  tree *tr = chain->tr; 
      
  proposalFunction
    *prevProposal = chain->prevProposal; 

  topoRecord
    *rec = thisProposal->remembrance.topoRec; 
  
  /* check, if the move changed anything at all */
  if( comparePair(rec->pruningBranches[0], rec->pruningBranches[1], rec->insertBranch[0], rec->insertBranch[1]) ) 
    {
      /* TODO here we could save tim e */      
      evaluateGenericWrapper(chain, chain->tr->start, TRUE); /* TODO   */      
      return ; 
    }

  /* cleanup necessary  */
  if(NOT chain->wasAccepted)
    {
      switch(prevProposal->category )
      	{
      	/* case TOPOLOGY: */
      	  /* evaluateGenericWrapper(chain, chain->tr->start, TRUE); */
      	  /* break; */
      	/* case RATE_HETEROGENEITY: */
      	  /* int modelChanged = prevProposal->remembrance.partInfo->modelNum;  */
      	  /* evaluatePartitions(chain, chain->tr->start, TRUE, ); */

      	default:
	  evaluateGenericWrapper(chain, chain->tr->start, TRUE); /* TODO   */
      	}
    }
  else 
    {       
      /* evaluate subtree again, if necessary   */
      nodeptr subtreePtr = getThirdNode(chain->tr, rec->prunedSubTree, rec->insertBranch[0], rec->insertBranch[1]);      
      exa_newViewGeneric( chain, subtreePtr, FALSE);

      /* find the ptr to the sub-tree that was left unchanged, when we pruned the subtree  */
      /* if(NOT isTip(rec->pruningBranches[1], tr->mxtips)) */
      /* 	{ */
      nodeptr goodSubtreePtr = findNodePtrByNeighbor( tr, rec->pruningBranches[1], rec->pruningBranches[0] );
      exa_newViewGeneric(chain, goodSubtreePtr, FALSE);
      /* 	} */

      /* evaluateGenericWrapper(chain, thisProposal->remembrance.topoRec->) */
      evaluateGenericWrapper(chain, chain->tr->start, TRUE); /* TODO   */
    }
#endif
}



static nodeptr select_random_subtree(state *chain, tree *tr)
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


static void extended_spr_apply(state *chain, proposalFunction *pf)
{
  tree *tr = chain->tr;
  topoRecord *rec = pf->remembrance.topoRec; 
  double stopProp = pf->parameters.eSprStopProb; 

  debug_printTree(chain);

  nodeptr prunedSubtreePtr = select_random_subtree(chain,tr);
  nodeptr nb = prunedSubtreePtr->next->back, 
    nnb = prunedSubtreePtr->next->next->back; 
  
  /* TODO can we assert this?  */
  assert( NOT isTip(prunedSubtreePtr->number, tr->mxtips) ) ; 

  rec->pruned = prunedSubtreePtr->number; 

  rec->pruningBranch.thisNode = nb->number; 
  rec->pruningBranch.thatNode = nnb->number; 

  record_branch_info(nb, rec->neighborBls , getNumBranches(chain->tr));
  record_branch_info(nnb, rec->nextNeighborBls, getNumBranches(chain->tr));
  
  double* nbz = rec->neighborBls,
    *nnbz = rec->nextNeighborBls; 

  /* initial remapping of BL of nodes adjacent to pruned node  */
  double zqr[NUM_BRANCHES];

  for(int i = 0; i < getNumBranches(tr); i++)
    {
      switch(pf->ptype)
	{
	case E_SPR : 
	  zqr [i] = nbz[i] * nnbz[i]; 
	  chain->hastings *= log(zqr[i]);
	  zqr[i] = sqrt(zqr[i]);
	  break; 
	case E_SPR_MAPPED: 
	  zqr[i] = nbz[i] ; 
	  break; 
	default: assert(0); 
	}

      if(zqr[i] > zmax) zqr[i] = zmax;
      if(zqr[i] < zmin) zqr[i] = zmin;
    }
  
  

  pruneBranch(chain,  constructBranch(prunedSubtreePtr->number, prunedSubtreePtr->back->number),zqr);
  /* pruneNodeFromNodes(chain, rec->pruned, rec->pruningBranch.thisNode, rec->pruningBranch.thatNode, zqr); */

  debug_printNodeEnvironment(chain, nb->number); 
  debug_printNodeEnvironment(chain, nnb->number); 

  /* done remove node p (omitted BL opt) */

  nodeptr initNode = NULL; 
  nodeptr curNode = NULL; 

  boolean remapBL = FALSE; 
  if(drawRandDouble01(chain) > 0.5)
    {
      curNode = nb; 
      remapBL = TRUE ; 
    }
  else 
    curNode = nnb; 
  initNode = curNode; 

  /* while(curNode->number == initNode->number */
  /* 	|| curNode->back->number == initNode->number) */
  /*   {				 */
      /* printf("initnode->numebr = %d\n", initNode->number);  */
  int accepted = FALSE;   
  while( NOT  accepted)
    {       
      accepted = extended_spr_traverse(chain, &curNode, stopProp ); 
      
      /* needed for spr remap */
      if(curNode == initNode)
	remapBL = NOT remapBL; 
    }
    /* } */

  rec->insertBranch.thisNode = curNode->number;   
  rec->insertBranch.thatNode = curNode->back->number; 


  nodeptr insertBranchPtr = curNode; 

  record_branch_info(insertBranchPtr, pf->remembrance.topoRec->bls, getNumBranches(chain->tr)); 

  switch(pf->ptype)
    {
      /* TODO this was the standard!  */
      /* case STANDARD:  */
      /* for(int branchCount=0; branchCount< getNumBranches(chain->tr); branchCount++) */
      /*   { */
      /*     chain->hastings/=log(curNode->z[branchCount]); */
      /*   }   */
      /*     insertWithUnifBL(chain->sprMoveRemem.p, chain->sprMoveRemem.q, getNumBranches(chain->tr)); */
      /* /\*   for(int branchCount=0; branchCount<chain->tr->numBranches; branchCount++) */
      /*    { */
      /*      chain->hastings/=log(chain->brLenRemem.qz[branchCount]); */
      /*    } */
      /*    *\/ */
      /*  //  printf("hastings: %f\n", chain->hastings); */
      /*     // insertBIG(chain->tr, chain->sprMoveRemem.p, chain->sprMoveRemem.q, chain->tr->numBranches); */
      break; 


    case E_SPR  :  		/* ADJUCT */
      for(int branchCount=0; branchCount< getNumBranches(chain->tr); branchCount++) /*  */
	{
	  chain->hastings/=(2*log(curNode->z[branchCount]));
	  insertWithUnifBLScaled(prunedSubtreePtr, insertBranchPtr, 2.0,  getNumBranches(chain->tr));
	}
      break;       
    case E_SPR_MAPPED: 
      {
	/* TODO hastings?  */
      
	double *neighborZ = nnbz; 

	/* if we went into the other direction, correct for that. 
	   specfically rec->pruningBranches must point into the direction, we moved the tree to 
	 */
	if( remapBL ) 
	  {
	    neighborZ = nbz; 
	    for(int i = 0; i < getNumBranches(tr); ++i)
	      nb->z[i] = nb->back->z[i] = nnb->z[i]; 
	    
	    assert(getNumBranches(tr) == 1 ); 
	    swpDouble(rec->neighborBls,rec->nextNeighborBls ); 
	    rec->pruningBranch = invertBranch(rec->pruningBranch); 
	  }
	

	insertNodeIntoBranch( chain, 
			      constructBranch(prunedSubtreePtr->number, prunedSubtreePtr->back->number ),
			      rec->insertBranch, 
			      insertBranchPtr->z,
			      neighborZ); 
	/* insertNodeIntoBranch(chain, prunedSubtreePtr->number, rec->insertBranch.thisNode,rec->insertBranch.thatNode, insertBranchPtr->z, neighborZ); */
      }
      break; 
    default: assert(0); 
    }

  debug_checkTreeConsistency(chain); 

  assert(branchExists(tr, constructBranch(rec->pruned, rec->insertBranch.thisNode ))); 
  assert(branchExists(tr, constructBranch(rec->pruned, rec->insertBranch.thatNode ))); 
} 



static void extended_spr_reset(state * chain, proposalFunction *pf)
{
  tree *tr = chain->tr; 
  topoRecord
    *topoRec = pf->remembrance.topoRec; 
  

  /* HACK to find the other node in the subtree */  
  nodeptr p = tr->nodep[topoRec->pruned] ; 
  while(p->back->number == topoRec->insertBranch.thisNode 
	|| p->back->number == topoRec->insertBranch.thatNode )
    p = p->next; 
  p = p->back; 
  assert(p->number != topoRec->insertBranch.thisNode  && p->number != topoRec->insertBranch.thatNode); 
  /* END */

  pruneBranch(chain, constructBranch(topoRec->pruned, p->number) , topoRec->bls);

  insertNodeIntoBranch(chain, 
		       constructBranch(topoRec->pruned, p->number), 
		       topoRec->pruningBranch, 
		       topoRec->neighborBls, topoRec->nextNeighborBls) ; 

#ifdef DEBUG_SHOW_TOPO_CHANGES
  debug_printNodeEnvironment(chain, topoRec->pruned ); 
  debug_printNodeEnvironment(chain,   topoRec->pruningBranches[0] ); 
  debug_printNodeEnvironment(chain,   topoRec->pruningBranches[1] ); 
  debug_printNodeEnvironment(chain,   topoRec->insertBranch[0] ); 
  debug_printNodeEnvironment(chain,   topoRec->insertBranch[1] ); 
#endif

  debug_checkTreeConsistency(chain );
  debug_printTree(chain);


  /* TODO this legacy stuff does not really make sense...  */
  /* exa_newViewGeneric(chain, prunedPtr, FALSE);  */
  /* double val1 = chain->tr->likelihood;  */  
  /* exa_newViewGeneric(chain, chain->tr->start, TRUE); */
  /* double  val2 = chain->tr->likelihood;  */
  /* assert( fabs ( val2 - val1 ) < 0.0001 );  */
}


//--------Alpha-Proposal-for-GAMMA-----------------------------------------------------------------

double get_alpha_prior(state *chain )
{
  
 return 1;//TODO obviously needs acctual prior 
}





static void simple_gamma_proposal_apply(state * chain, proposalFunction *pf)
{
  tree *tr = chain->tr; 

  //TODO: add safety to max and min values
  double newalpha, curv, r,mx,mn;

  int model = drawRandInt(chain, getNumberOfPartitions(tr));  
  perPartitionInfo* remem = pf->remembrance.partInfo; 
  remem->modelNum = model; 

  pInfo *partition = getPartition(chain,model);
  curv = partition->alpha;
  remem->alpha = curv; 

  switch(pf->ptype)
    {
    case UPDATE_GAMMA:
      {
	double slidWin = pf->parameters.slidWinSize;
	    
	/* case STANDARD://simple sliding window */
	r = drawRandDouble01(chain);
	mn = curv-(slidWin/2);
	mx = curv+(slidWin/ 2);
	newalpha = fabs(mn + r * (mx-mn));
      }
      break;
    case UPDATE_GAMMA_EXP: 
      newalpha  = drawRandExp(chain,1/curv);
      break;
    default:
      assert(0);
    }
  /* Ensure always you stay within this range */
  if(newalpha > ALPHA_MAX) newalpha = ALPHA_MAX;
  if(newalpha < ALPHA_MIN) newalpha = ALPHA_MIN;
  
  chain->hastings = 1; //since it is symmetrical, hastings=1
  
  partition->alpha = newalpha;
  makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);

  /* evaluateOnePartition(chain, tr->start, TRUE, model );  */
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


static void simple_gamma_proposal_reset(state *chain, proposalFunction *pf)
{
  tree *tr = chain->tr; 
  perPartitionInfo *info  = pf->remembrance.partInfo; 
  pInfo *partition = getPartition(chain, info->modelNum) ; 
  
  partition->alpha = info->alpha; 

  makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);
 /* evaluateOnePartition(chain, tr->start, TRUE, info->modelNum);  */
}

//------------------------------------------------------------------------------


#if 0 
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

#endif





//NOTE: should only be called at the very beginning. Afterwards a probability sum of 1.0 is not required.


#if  0
void normalizeProposalWeights(state *chain)
{
  double sum = 0 ; 
  for(int i = 0; i < NUM_PROPOSALS;++i)
    sum += chain->proposalWeights[i]; 
  
  for(int i = 0; i < NUM_PROPOSALS;++i)
    chain->proposalWeights[i] /= sum; 
  
  findLogisticT(chain);
}
#endif



static void simple_model_proposal_apply(state *chain, proposalFunction *pf)//llpqr
{

  //TODO: add safety to max and min values
  //record the old ones

  int model = drawRandInt(chain, getNumberOfPartitions(chain->tr));
  perPartitionInfo *info = pf->remembrance.partInfo;

  pInfo *partition = getPartition(chain, model); 
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
  
  int list[numRates];//for biunif_distr and biunif_perm_distr
  
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
	    curv = partition->substRates[state];
	    r =  drawRandDouble01(chain);
	    mn = curv-range ;
	    mx = curv+range ;
	
	    new_value = fabs(mn + r * (mx-mn));
	  
	    /* Ensure always you stay within this range */
	    if(new_value > RATE_MAX) new_value = RATE_MAX;
	    if(new_value < RATE_MIN) new_value = RATE_MIN;
	  }
	  break;
        
	case UPDATE_MODEL_BIUNIF:
	  changeState=drawRandInt(chain, numRates);
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

  exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
    
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



static void perm_biunif_model_proposal_apply(state *chain, proposalFunction *pf)
{
  /* tree *tr = chain->tr;  */

  int model = drawRandInt(chain,getNumberOfPartitions(chain->tr));; 

  pInfo *partition = getPartition(chain , model) ; 

  perPartitionInfo *info = pf->remembrance.partInfo; 

  int numRates  = (partition->states  * partition->states -partition->states ) / 2 ;  
  info->numRates = numRates; 

  recordSubsRates(chain, model, numRates, info->substRates);
  int state, randNumber;
  double new_value,curv;
  double r;
  
  randNumber=drawRandInt(chain,numRates);
  int perm[numRates];
  drawPermutation(chain,perm, numRates);

  for(state = 0;state < randNumber ; state ++)
    {           
      curv = partition->substRates[perm[state]];
      r =  drawRandBiUnif(chain,curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
	if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
	if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
      chain->hastings*=curv/new_value;
      edit_subs_rates(chain,model, perm[state], new_value);      
    }
      
  exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/   */
}



static void single_biunif_model_proposal_apply(state *chain,proposalFunction *pf)//NOTE whenever a model parameter changes, all branch lengths have to be re-normalized with 1/fracchange. Additionally we always must do a full tree traversal to get the likelihood. So updating a single parameter is rather expensive, .
{
  /* tree *tr = chain->tr;  */

  //record the old one //TODO sufficient to store single value.  
  int model = drawRandInt(chain,getNumberOfPartitions(chain->tr)); 

  perPartitionInfo *info =  pf->remembrance.partInfo; 
  
  pInfo *partition = getPartition(chain,model) ; 

  info->modelNum = model; 
  info->numRates = partition->states ;   
  int numRates = info->numRates; 
  recordSubsRates(chain, model, info->numRates, info->substRates);
  //choose a random set parameter,
  //with uniform probabilities

  int  randState=drawRandInt(chain,numRates);

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

  edit_subs_rates(chain,model, randState, new_value);

  chain->hastings=curv/new_value;

  exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
  
  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/ */
}

static void all_biunif_model_proposal_apply(state *chain, proposalFunction *pf)
{
  /* tree *tr = chain->tr;  */
  
  int model = drawRandInt(chain,getNumberOfPartitions(chain->tr));

  perPartitionInfo *info = (perPartitionInfo* ) pf->remembrance.partInfo; 
  
  info->modelNum = model; 
  
  //record the old one 
  pInfo *partition = getPartition(chain, model) ; 
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
      curv = partition->substRates[state]; 
      r =  drawRandBiUnif(chain,curv);

      new_value = r;

      while(new_value> RATE_MAX|| new_value< RATE_MIN){
      if(new_value > RATE_MAX) new_value = 2*RATE_MAX-new_value;
      if(new_value< RATE_MIN) new_value= 2*RATE_MIN-new_value;
      }
      
       chain->hastings*=curv/new_value;
      edit_subs_rates(chain,model, state, new_value);
    }

  exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); /\* 2. re-traverse the full tree to update all vectors *\/ */
}

static void restore_subs_rates(state *chain, int model, int numSubsRates, double *prevSubsRates)
{
  /* tree *tr = chain->tr;    */

  pInfo *partition = getPartition(chain, model); 

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)	
    partition->substRates[i] = prevSubsRates[i]; 

  exa_initReversibleGTR(chain, model);

  /* TODO need to broadcast rates here for parallel version */

  /* evaluateOnePartition(chain, tr->start, TRUE, model);  */
}


//--------Branch-Length_Proposals---------------------------------------------------

double get_branch_length_prior( state *chain)
{//TODO decide on sensible prior
  return 1;  
}

//setting this out to allow for other types of setting
static void set_branch_length_sliding_window(state *chain, nodeptr p, int numBranches,state * s, boolean record_tmp_bl, double windowRange, double *bls)
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
	  bls[i] = p->z[i]; 

	  /* bls[i] = p->back->z_tmp[i] = p->z[i]; */
	  /* p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /\* keep current value *\/ */
	}
      r = drawRandDouble01(chain);

      /* TODO@kassian just wondering: is it correct to use the global
	 frac-change here? what are the local frac-changes good for
	 then?  */
      real_z = -log(p->z[i]) * s->tr->fracchange;
      
      //     printf( "z: %f %f\n", p->z[i], real_z );
    
      mn = real_z- windowRange;
      mx = real_z+ windowRange;
      newZValue = exp(-(fabs(mn + r * (mx-mn)/s->tr->fracchange )));

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

static void set_branch_length_biunif(state *chain, nodeptr p, int numBranches,state * s, boolean record_tmp_bl, double *bls)
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
	  bls[i] = p->z[i]; 
	  /* p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /\* keep current value *\/ */
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

static void set_branch_length_exp(state *chain, nodeptr p, int numBranches,state * s, boolean record_tmp_bl, double *bls)
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
	  bls[i] = p->z[i]; 
	  /* p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /\* keep current value *\/ */
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




static void random_branch_length_proposal_apply(state * chain, proposalFunction *pf)
{
  int multiBranches = getNumBranches(chain->tr);
  int target_branch = drawRandInt(chain,(chain->tr->mxtips * 2) - 3); 
  node *p = select_branch_by_id_dfs( chain->tr->start, target_branch, chain );
  
  switch(pf->ptype)
    {
    case UPDATE_SINGLE_BL: 
      //for each branch get the current branch length
      //pull a uniform like
      //x = current, w =window
      //uniform(x-w/2,x+w/2)
      set_branch_length_sliding_window(chain,p, multiBranches, chain, TRUE, pf->parameters.slidWinSize, pf->remembrance.topoRec->bls); 
      break;
	
    case UPDATE_SINGLE_BL_EXP: 
      set_branch_length_exp(chain,p, multiBranches, chain, TRUE,pf->remembrance.topoRec->bls);
      break;
    case UPDATE_SINGLE_BL_BIUNIF: 

      set_branch_length_biunif(chain, p, multiBranches, chain, TRUE,pf->remembrance.topoRec->bls);
      break;
    default:
      assert(0);
    }

  pf->remembrance.topoRec->whichBranch = target_branch; 
  /* evaluateGenericWrapper(chain, p, FALSE);  */
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

static void random_branch_length_proposal_reset(state * chain, proposalFunction *pf)
{
  node *p;

  // ok, maybe it would be smarter to store the node ptr for rollback rather than re-search it...
  p = select_branch_by_id_dfs( chain->tr->start, pf->remembrance.topoRec->whichBranch, chain );
  
  reset_branch_length(p, getNumBranches(chain->tr), pf->remembrance.topoRec->bls);
  //   printf( "reset bl: %p %f\n", p, p->z[0] );
  //update_all_branches(chain, TRUE);

  /* i am not so sure, if we really can do the FALSE here ; will raxml recognize that a branch length has changed? disabling it for now... 

     TODO I think, we should evaluate at the respctive node 
   */

 // evaluateGeneric(chain->tr, p, FALSE); //This yields a very slight likelihood difference.NOTE if we want exact likelihoods as before the proposal, we must evaluate from chain->tr->start, that is: evaluateGeneric(chain->tr, chain->tr->start, TRUE);
}

double get_frequency_prior(state * chain)
{
 return 1; 
}

static void restore_frequ_rates(state *chain, int model, int numFrequRates, double *prevFrequRates)
{
  /* tree *tr = chain->tr; */

  /* NOTICE: this function should not be called repeatedly  */

  pInfo *partition = getPartition(chain,model);

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    partition->frequencies[i] = prevFrequRates[i];

  exa_initReversibleGTR(chain, model);

  /* evaluateOnePartition(chain, tr->start, TRUE, model); */
}

static void recordFrequRates(state *chain, int model, int numFrequRates, double *prevFrequRates)
{
  pInfo *partition = getPartition(chain,model);

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    prevFrequRates[i] = partition->frequencies[i];
}




void frequency_proposal_apply(state * chain, proposalFunction *pf)
{
  tree *tr = chain->tr; 

  int model  = drawRandInt(chain,getNumberOfPartitions(tr));
  perPartitionInfo *info = pf->remembrance.partInfo; 

  pInfo *partition = getPartition(chain, model); 

  int numFreq = partition->states; 
  info->modelNum = model;   
  info->numFreqs = numFreq; 

  recordFrequRates(chain, model, numFreq, info->frequencies);
  
  double r[numFreq];  
  for(int state = 0;state < numFreq ; state ++)
    {
      double curv = partition->frequencies[state];
      r[state] = drawRandBiUnif(chain,curv);       
      chain->hastings*=curv/r[state];      
    }

  double sum=0;  
  for(int state = 0;state< numFreq ; state ++)    
    sum+=r[state]; 

  for(int state = 0; state< numFreq ; state ++)    
    partition->frequencies[state]=r[state]/sum; 

  //recalculate eigens

  exa_initReversibleGTR(chain, model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* chain->curprior=get_frequency_prior(chain, tr->partitionData[chain->frequRemem.model].frequencies); */
  /* chain->newprior=get_frequency_prior(chain, chain->frequRemem.curFrequRates); */

  /* evaluateOnePartition(chain, tr->start, TRUE, model); */
}


void frequency_proposal_reset(state * chain, proposalFunction *pf)
{
  perPartitionInfo *info = pf->remembrance.partInfo; 
  restore_frequ_rates(chain, info->modelNum, info->numFreqs, info->frequencies);
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



static void simple_model_proposal_reset(state * chain, proposalFunction *pf)
{
  perPartitionInfo *info = pf->remembrance.partInfo; 
  restore_subs_rates(chain, info->modelNum, info->numRates, info->substRates);
}



static void branch_length_multiplier_apply(state *chain, proposalFunction *pf)
{
  tree *tr = chain->tr; 
  branch b =  drawBranchUniform(chain); 

  nodeptr p = findNodeFromBranch(tr, b); 

  topoRecord
    *rec = pf->remembrance.topoRec; 
  rec->insertBranch = b ; 

  record_branch_info( p, rec->bls, getNumBranches(tr));

  double
    multiplier = exp(pf->parameters.multiplier * (drawRandDouble01(chain) - 0.5 ));   
  assert(multiplier > 0.); 
  
  /* TODO how do we do that wiht multiple bls per branch?  */
  assert(getNumBranches(tr) == 1); 
  
  double newZ  = multiplier * ( -log(p->z[0]) * tr->fracchange);  
  double newZTransformed = exp(-(newZ / tr->fracchange)); 
  
  /* according to lakner2008  */
  chain->hastings *= multiplier; 
 
  /* just doing it for one right here */
  p->z[0] = p->back->z[0] = newZTransformed; 
}


/**
   @brief tunes the BL multiplier
 */ 
static void autotuneMultiplier(state *chain, proposalFunction *pf)
{
  double *parameter = &(pf->parameters.multiplier); 

  successCtr *ctr = &(pf->sCtr); 

  int batch = chain->currentGeneration  / TUNE_FREQUENCY; /* HACK */

  double newParam = tuneParameter(batch, getRatioLocal(ctr), *parameter, FALSE); 

#ifdef DEBUG_PRINT_TUNE_INFO
  if(processID == 0 )
    {
    if(newParam < *parameter)
      printf("%s\tratio=%f\t => reduced %f to %f\n", pf->name, getRatioLocal(ctr), *parameter, newParam);
    else if (newParam > *parameter)
      printf("%s\tratio=%f\t => increased %f to %f\n", pf->name, getRatioLocal(ctr), *parameter, newParam);
    }
#endif

  *parameter = newParam; 
  resetCtr(ctr);   
}


/**
   @brief autotunes sliding windows 
   
*/ 
static void autotuneSlidingWindow(state *chain, proposalFunction *pf)
{
  double *parameter = &(pf->parameters.slidWinSize); 
  successCtr *ctr = &(pf->sCtr); 
  double newParam = tuneParameter(chain->currentGeneration / TUNE_FREQUENCY,
				  getRatioLocal(&(pf->sCtr)),
				  *parameter, FALSE  ); 
  
#ifdef DEBUG_PRINT_TUNE_INFO
  if(processID == 0 )
    {
      if(newParam < *parameter)
	printf("%s ratio=%f\t => reduced %f to %f\n", pf->name, getRatioLocal(ctr), *parameter, newParam);
      else if (newParam > *parameter)
	printf("%s ratio=%f\t => increased %f to %f\n", pf->name, getRatioLocal(ctr), *parameter, newParam);
    }
#endif

  *parameter = newParam; 
  resetCtr(ctr);   

  

}


/**
   @brief autotunes the stop probability of extended topological moves 
*/ 
static void autotuneStopProp(state *chain, proposalFunction *pf) 
{
  double *parameter = &(pf->parameters.eSprStopProb); 
  successCtr *ctr = &(pf->sCtr); 
  double newParam = tuneParameter( chain->currentGeneration / TUNE_FREQUENCY, 
				   getRatioLocal(&(pf->sCtr)),
				   *parameter, 
				   TRUE ); 
  
#ifdef DEBUG_PRINT_TUNE_INFO
  if(processID == 0 )
    {
    if(newParam < *parameter)
      printf("%s ratio=%f\t => reduced %f to %f\n", pf->name,getRatioLocal(ctr), *parameter, newParam);
    else if (newParam > *parameter)
      printf("%s ratio=%f\t => increased %f to %f\n", pf->name, getRatioLocal(ctr), *parameter, newParam);
    }
#endif

  *parameter = fmin(newParam,1); 
  resetCtr(ctr);   
  
}


void branch_length_reset(state *chain, proposalFunction *pf)
{
  tree *tr = chain->tr; 
  topoRecord *rec = pf->remembrance.topoRec; 
  nodeptr p = findNodeFromBranch(tr, rec->insertBranch); 
  
  reset_branch_length(p,getNumBranches(tr), rec->bls); 

  /* evaluateGenericWrapper(chain, tr->start, TRUE );    */
}



static void initProposalFunction( proposal_type type, initParamStruct *initParams, proposalFunction **result)
{
  if(initParams->initWeights[type] == 0)
    {
      *result = NULL; 
      return ; 
    }  

  *result = exa_calloc(1,sizeof(proposalFunction));   

  proposalFunction *ptr = *result; 
  ptr->ptype = (proposal_type)type; 
  ptr->initWeight = initParams->initWeights[type]; 
  ptr->currentWeight = ptr->initWeight; 


  /* 
     TODO@kassian
   
     I fear I could not restore all proposals correctly. Or let me put
     it differently, I do not trust all of them yet. Those, that I do
     trust are marked.
  */

  switch(type)
    {
    case E_SPR:
      ptr->eval_lnl = dummy_eval;
      ptr->autotune = autotuneStopProp; 
      ptr->remembrance.topoRec = exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func = extended_spr_apply; 
      ptr->reset_func = extended_spr_reset; 
      ptr->category = TOPOLOGY;       
      ptr->parameters.eSprStopProb = initParams->eSprStopProb; 
      ptr->name = "eSPR"; 
      break; 
    case E_SPR_MAPPED: 		/* TRUSTED  */
      ptr->eval_lnl = dummy_eval;
      ptr->autotune = autotuneStopProp; 
      ptr->remembrance.topoRec = exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func = extended_spr_apply; 
      ptr->reset_func = extended_spr_reset; 
      ptr->name = "eSPRMapped"; 
      ptr->parameters.eSprStopProb = initParams->eSprStopProb; 
      ptr->category = TOPOLOGY; 
      break; 
    case UPDATE_MODEL: 	
      ptr->eval_lnl = onePartitionEval; 
      /* ptr->eval_lnl = dummy_eval; */
      ptr->autotune = autotuneSlidingWindow; 
      ptr->apply_func = simple_model_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->parameters.slidWinSize = INIT_RATE_SLID_WIN;
      ptr->category = SUBSTITUTION_RATES; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "modelSlidWin"; 
      break; 
    case UPDATE_GAMMA:      	
      ptr->eval_lnl = onePartitionEval;
      /* ptr->eval_lnl = dummy_eval; */
      ptr->autotune = autotuneSlidingWindow; 
      ptr->apply_func = simple_gamma_proposal_apply; 
      ptr->reset_func = simple_gamma_proposal_reset; 
      ptr->parameters.slidWinSize = INIT_RATE_SLID_WIN; 
      ptr->category = RATE_HETEROGENEITY; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "gammaSlidWin"; 
      break; 
    case UPDATE_GAMMA_EXP: 
      /* ptr->eval_lnl = dummy_eval; */
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = simple_gamma_proposal_apply; 
      ptr->reset_func = simple_gamma_proposal_reset; 
      ptr->category = RATE_HETEROGENEITY; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "gammaExp"; 
      break; 
    case UPDATE_SINGLE_BL: 	/* TRUSTED */
      ptr->eval_lnl = dummy_eval;
      ptr->autotune = autotuneSlidingWindow; 
      ptr->remembrance.topoRec = exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func	=  random_branch_length_proposal_apply;
      ptr->reset_func =  random_branch_length_proposal_reset;
      ptr->category = BRANCH_LENGTHS; 
      ptr->parameters.slidWinSize = INIT_BL_SLID_WIN; 
      ptr->name = "singleBLSlidWin"; 
      break; 
    case UPDATE_SINGLE_BL_EXP: 
      ptr->eval_lnl = dummy_eval;
      ptr->remembrance.topoRec = exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func	=  random_branch_length_proposal_apply;
      ptr->reset_func =  random_branch_length_proposal_reset;
      ptr->category = BRANCH_LENGTHS; 
      ptr->name  = "singleBlExp"; 
      break; 
    case UPDATE_SINGLE_BL_BIUNIF: 
      ptr->eval_lnl = dummy_eval;
      ptr->remembrance.topoRec = exa_calloc(1,sizeof(topoRecord)); 
      ptr->apply_func	=  random_branch_length_proposal_apply;
      ptr->reset_func =  random_branch_length_proposal_reset;
      ptr->name = "singleBLBiunif"; 
      ptr->category = BRANCH_LENGTHS; 
      break; 
    case UPDATE_MODEL_SINGLE_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      /* ptr->eval_lnl = dummy_eval; */
      ptr->apply_func	=  single_biunif_model_proposal_apply;
      ptr->reset_func =  simple_model_proposal_reset;
      ptr->name = "singleModelBiunif"; 
      ptr->category = SUBSTITUTION_RATES; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      break; 
    case UPDATE_MODEL_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      /* ptr->eval_lnl = dummy_eval; */
      ptr->name = "modelBiunif"; 
      ptr->category = SUBSTITUTION_RATES; 
      ptr->apply_func = single_biunif_model_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      break; 
    case UPDATE_MODEL_ALL_BIUNIF:
      ptr->eval_lnl = onePartitionEval; 
      /* ptr->eval_lnl = dummy_eval;  */
      ptr->apply_func = all_biunif_model_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->name = "modelAllBiunif"; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->category = SUBSTITUTION_RATES; 
      break; 
    case UPDATE_FREQUENCIES_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      /* ptr->eval_lnl = dummy_eval; */
      ptr->apply_func = frequency_proposal_apply; 
      ptr->reset_func = frequency_proposal_reset; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "freqBiunif"; 
      ptr->category = FREQUENCIES; 
      break;
    case UPDATE_MODEL_PERM_BIUNIF: 
      ptr->eval_lnl = onePartitionEval; 
      /* ptr->eval_lnl = dummy_eval; */
      ptr->apply_func = NULL; 	/* TODO */
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->remembrance.partInfo = exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "modelPermBiunif"; 
      ptr->category = SUBSTITUTION_RATES; 
      break;
    case BRANCH_LENGTHS_MULTIPLIER: 
      ptr->eval_lnl = dummy_eval;
      ptr->autotune = autotuneMultiplier; 
      ptr->apply_func = branch_length_multiplier_apply; 
      ptr->reset_func = branch_length_reset; 
      ptr->remembrance.topoRec = exa_calloc(1,sizeof(topoRecord)); 
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






/**
   @brief Normalizes the weights of the proposals in this category
 */
void normalizePropSubCats(state *chain)
{
  double* catWeights = exa_calloc(NUM_PROP_CATS + 1,sizeof(double)); 

  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      assert(pf->category); 
      catWeights[pf->category] +=  pf->currentWeight;       
    }

  for(int i= 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      pf->currentWeight /= catWeights[pf->category]; 
    }

  exa_free(catWeights); 
}


/**
   @brief normalizes the categories  
 */
void normalizeCategories(state *chain)
{
  double sum = 0.; 
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    sum += chain->categoryWeights[i]; 
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    chain->categoryWeights[i] /= sum;   
}




void printAllProposalWeights(state *chain)
{
  if(processID != 0)
    return; 
  
  printf("cat weights: TOPO=%f\tBL=%f\tFREQ=%f\tSUBST=%f\tHET=%f\n", 
	 chain->categoryWeights[0],
	 chain->categoryWeights[1],
	 chain->categoryWeights[2],
	 chain->categoryWeights[3],
	 chain->categoryWeights[4] ); 

  printf("rel. prop weihgts: "); 
  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      printf("\t%s=%f", pf->name, pf->currentWeight);       
    }
  printf("\n"); 
}



/**
   @brief   Initializes the proposals based on weights given in the config file. 

   Also normalizes all weights to 1.   

 */
void setupProposals(state *chain, initParamStruct *initParams)
{
  int ctr = 0; 

  proposalFunction **pfs = exa_calloc(NUM_PROPOSALS, sizeof(proposalFunction*));  
  for(int i = 0; i < NUM_PROPOSALS; ++i)
    {
      proposalFunction *pf = NULL; 
      initProposalFunction((proposal_type)i, initParams, &pf); 
      if(pf != (proposalFunction*)NULL)
	pfs[ctr++] = pf; 
    }  
  chain->numProposals = ctr; 
  chain->proposals = pfs; 

  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      chain->categoryWeights[pf->category-1] += pf->currentWeight; 
    }
  normalizeCategories(chain);  
  normalizePropSubCats(chain); 
  
  /* only print that once  */
  if(chain->id == 0)
    printAllProposalWeights(chain);
}



/**
   @brief draws a proposal function.

   Notice: this could be extended later, if we decide to make this
   dependent on the previous state.
   
   Furthermore, we must be sure now that category weights and relative
   proposal weights sum up to 1 each. 
   
 */ 
void drawProposalFunction(state *chain, proposalFunction **result )
{  
  
  *result = NULL; 
  category_t
    cat = drawSampleProportionally(chain,chain->categoryWeights, NUM_PROP_CATS) + 1; /* it is 1-based */

  /* printInfo(chain, "drawing proposal; category is %d\n"), cat;  */
  
  double sum = 0; 
  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      if(pf->category == cat)
	sum += pf->currentWeight; 
    }
  assert(fabs(sum - 1.) < 0.000001); 

  double r = drawRandDouble01(chain);
  /* printf("numProp=%d\n", chain->numProposals);  */
  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      if(pf->category == cat)
	{
	  if(  r < pf->currentWeight)
	    {
	      *result =  pf; 
	      return; 
	    }
	  else 
	    {
	      r -= pf->currentWeight; 
	    }
	}
    }

  assert(result != NULL); 
}




void debug_checkLikelihood(state *chain)
{
#ifdef DEBUG_LNL_VERIFY
  tree *tr = chain->tr; 
  exa_evaluateGeneric(chain, tr->start, TRUE );   
  double diff =  fabs(chain->likelihood - tr->likelihood ); 
  if( diff > ACCEPTED_LIKELIHOOD_EPS)
    {
      printf("LNL difference: %f\n", diff  ); 
      assert(diff  < ACCEPTED_LIKELIHOOD_EPS); 
    }  
  chain->likelihood = tr->likelihood; 
#endif  
}







/**
   @brief Execute one generation of a given chain.  
 */
void step(state *chain)
{
  tree *tr = chain->tr;   

  double prevLnl = chain->lnl.likelihood;    
  int numPart = getNumberOfPartitions( chain->tr );
  double prevPartLnl[numPart]; 
  for(int i = 0; i < numPart; ++i)
    prevPartLnl[i] = chain->lnl.partitionLnl[i]; 

  double myHeat = getChainHeat(chain ) ; 

  proposalFunction *pf = NULL;   
  drawProposalFunction(chain, &pf);

  /* reset proposal ratio  */
  chain->hastings = 1; 

  double oldPrior = chain->priorProb; 		/* TODO  */

  /* chooses move, sets proposal ratio, correctly modifies the prior */
  pf->apply_func(chain, pf);  
  double priorRatio  = chain->priorProb - oldPrior; 
  /* enable once we actually have priors  */
  assert(priorRatio == 0); 

  /* chooses the cheapest way to evaluate the likelihood  */
  pf->eval_lnl(chain, pf); 

  double testr = drawRandDouble01(chain);
  double acceptance = 
    exp((priorRatio  + chain->lnl.likelihood - prevLnl) * myHeat) 
    * chain->hastings ; 

  chain->wasAccepted  = testr < acceptance; 
  debug_printAccRejc(chain, pf, chain->wasAccepted); 
  chain->prevProposal = pf;   

  if(chain->wasAccepted)
    {
      cntAccept(&(pf->sCtr));
      
      double tmp = 0; 
      for(int i = 0; i < numPart; ++i)
	tmp += chain->lnl.partitionLnl[i]; 
      assert(fabs(tmp - chain->lnl.likelihood) < 1e-6); 

      /* 
	 commenting this out for now because of drastic changes: but
	 we should re-enable that later and do a category-wide tuning
	 of moves	 
       */
      /* penalize(chain, which_proposal, 1); */
    }
  else
    {
      pf->reset_func(chain, pf); 
      cntReject(&(pf->sCtr)); 
      
      chain->lnl.likelihood = prevLnl; 
      double tmp = 0;
      for(int i = 0; i < numPart; ++i)
	{
	  chain->lnl.partitionLnl[i] = prevPartLnl[i]; 
	  tmp += prevPartLnl[i];
	}
      if(fabs(tmp - chain->lnl.likelihood) > 1e-6)
	{
	  printf("WARNING: lnl diff=%f\n", tmp - chain->lnl.likelihood); 
	  assert(0); 
	}
      
      /* TODO re-enable */
      /* penalize(chain, which_proposal, 0); */
    }

  debug_checkTreeConsistency(chain);
  debug_checkLikelihood(chain); 

  if(  processID == 0 
       && chain->couplingId == 0	/* must be the cold chain  */
       && (chain->currentGeneration % gAInfo.samplingFrequency) == gAInfo.samplingFrequency - 1  ) 
    {
      printSample(chain);       

      if(chain->currentGeneration >  BURNIN)
	addBipartitionsToHash(tr, chain ); 
    }

  
  /* the output for the console  */
  if(processID == 0 
     && chain->couplingId == 0     
     && chain->currentGeneration % PRINT_FREQUENCY == PRINT_FREQUENCY -1   )
    {
      chainInfo(chain); 
    }


#ifdef TUNE_PARAMETERS
  /* autotuning for proposal parameters. With increased parallelism
     this will become more complicated.  */
  if( chain->currentGeneration % TUNE_FREQUENCY == TUNE_FREQUENCY - 1 )
    {
      for(int i = 0; i < chain->numProposals; ++i)
	{
	  proposalFunction *pf = chain->proposals[i]; 
	  
	  if(pf->autotune)	/* only, if we set this   */
	    {
	      pf->autotune(chain, pf);
	    }
	}
    }
#endif

  chain->currentGeneration++; 
}

