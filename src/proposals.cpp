
#include "axml.h"
#include "bayes.h"		// 
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"
#include "topology-utils.h" 
#include "eval.h"
#include "adapters.h"
#include "misc-utils.h"
#include "branch.h"
#include "path.h"
#include "guidedMoves.h"
#include "stNNI.h"
#include "tlMult.h" 
#include "TreeAln.hpp"
#include "nodeSlider.h"

void expensiveVerify(tree *tr); 

/* nodeptr select_random_subtree(Chain *chain, tree *tr); */
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
  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = chain->traln->getSubstRate(model,i); 
}


// static void reset_branch_length(nodeptr p, int numBranches, double *bls)
// {
//   int i;
//   for(i = 0; i < numBranches; i++) 
//     p->z[i] = p->back->z[i] = bls[i];   /* restore saved value */
// }


static void record_branch_info(TreeAln* traln, nodeptr p, double *bl)
{
  int i;
  int numBranches = traln->getNumBranches(); 
  for(i = 0; i < numBranches; i++)
    bl[i] = traln->getBranchLength( p,0);
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
int extended_spr_traverse(Chain *chain, nodeptr *insertNode, double stopProp)
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


static void onePartitionEval(Chain *chain, proposalFunction *thisProposal)
{
  int model = thisProposal->remembrance.partInfo->modelNum;
  branch root = findRoot(chain->traln->getTr());  
  nodeptr p = findNodeFromBranch(chain->traln->getTr(), root); 
  evaluateOnePartition(chain, p, TRUE, model); 
}


static void evalBranch(Chain *chain, proposalFunction *thisProposal)
{
  branch b = peekStack(thisProposal->remembrance.modifiedPath);
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


void pushToStackIfNovel(stack *s, branch b, int numTip); 



static void sprEval(Chain *chain, proposalFunction *thisProposal)
{  
  tree *tr = chain->traln->getTr(); 
  path *rPath = thisProposal->remembrance.modifiedPath; 

  branch futureRoot = getThirdBranch(tr, rPath->content[0], rPath->content[1]); 
  
  /* evaluate at root of inserted subtree */
  nodeptr toEval = findNodeFromBranch(tr, futureRoot); /* dangerous */

  destroyOrientationAlongPath(tr, rPath, toEval); 
  destroyOrientationAlongPath(tr,rPath, toEval->back);

  /* printf("evaluating at branch %d,%d\n", toEval->number, toEval->back->number);  */
  evaluateGenericWrapper(chain, toEval, FALSE);
}



static void resetESPR(Chain *chain, proposalFunction *pf )
{
  resetAlongPathForESPR (chain->traln, pf->remembrance.modifiedPath);   

  /* TODO resetAlong... should be able to correctly restore branch lengths    */
  restoreBranchLengthsPath(chain->traln, pf->remembrance.modifiedPath); 

  debug_checkTreeConsistency(chain->traln->getTr());
  freeStack(&(pf->remembrance.modifiedPath));
  /* printf("RESET: ");    */
  debug_printTree(chain); 
}




/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
static void applyExtendedSPR(Chain *chain, proposalFunction *pf)
{
  double stopProp = pf->parameters.eSprStopProb; 

  debug_printTree(chain);
  
  path *rPath = NULL; 
  createStack(&rPath); 
  drawPathForESPR(chain,rPath,stopProp); 

  /* printStack(rPath); */

  saveBranchLengthsPath(chain, rPath); 

  applyPathAsESPR(chain->traln, rPath);

  pf->remembrance.modifiedPath = rPath; 

#ifdef ESPR_MULTIPLY_BL
  /* if(drawRandDouble01(chain) < 0.5) */
  multiplyAlongBranchESPR(chain, rPath, pf->param2.multiplier);
#endif

  debug_checkTreeConsistency(chain->traln->getTr()); 
}







static void extended_spr_apply(Chain *chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 
  tree *tr = chain->traln->getTr();
  topoRecord *rec = pf->remembrance.topoRec; 
  double stopProp = pf->parameters.eSprStopProb; 

  debug_printTree(chain);

  nodeptr prunedSubtreePtr = select_random_subtree(chain,tr);
  nodeptr nb = prunedSubtreePtr->next->back, 
    nnb = prunedSubtreePtr->next->next->back; 

  assert( NOT isTip(prunedSubtreePtr->number, tr->mxtips) ) ; 

  rec->pruned = constructBranch(0,prunedSubtreePtr->number); 
  rec->pruningBranch.thisNode = nb->number; 
  rec->pruningBranch.thatNode = nnb->number; 

  record_branch_info(traln, nb, rec->neighborBls);
  record_branch_info(traln,nnb, rec->nextNeighborBls);
  
  double* nbz = rec->neighborBls,
    *nnbz = rec->nextNeighborBls; 

  /* initial remapping of BL of nodes adjacent to pruned node  */
  double zqr[NUM_BRANCHES];

  for(int i = 0; i < chain->traln->getNumBranches(); i++)
    {
      switch(pf->ptype)
	{
	case E_SPR: 
	  zqr[i] = nbz[i] ; 
	  break; 
	default: assert(0); 
	}

      if(zqr[i] > zmax) zqr[i] = zmax;
      if(zqr[i] < zmin) zqr[i] = zmin;
    }
  
  pruneBranch(chain, constructBranch(prunedSubtreePtr->number, prunedSubtreePtr->back->number),zqr);

  debug_printNodeEnvironment(chain, nb->number); 
  debug_printNodeEnvironment(chain, nnb->number); 

  nodeptr initNode = NULL; 
  nodeptr curNode = NULL; 

  boolean remapBL = FALSE; 
  if(drawRandDouble01(chain) > 0.5)
    {
      curNode = nb; 
      remapBL = TRUE ; 
    }
  else 
    {
      curNode = nnb; 
    }
  initNode = curNode; 

  stack *s = NULL; 
  createStack(&s); 

  int accepted = FALSE;   
  while( NOT  accepted)
    {       
      accepted = extended_spr_traverse(chain, &curNode, stopProp ); 

      pushToStackIfNovel(s,constructBranch(curNode->number, curNode->back->number), tr->mxtips);

      /* needed for spr remap */
      if(curNode == initNode)
	{
	  remapBL = NOT remapBL; 
	  accepted = FALSE; 
	}

      if(stackIsEmpty(s))
	accepted = FALSE; 
    }
  
  /* printStack(s);  */
  freeStack(&s); 

  rec->insertBranch.thisNode = curNode->number;   
  rec->insertBranch.thatNode = curNode->back->number; 

  nodeptr insertBranchPtr = curNode; 

  record_branch_info(traln, insertBranchPtr, pf->remembrance.topoRec->bls); 

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

      break;       
    case E_SPR: 
      {
	/* TODO hastings?  */
      
	double *neighborZ = nnbz; 

	/* if we went into the other direction, correct for that. 
	   specfically rec->pruningBranches must point into the direction, we moved the tree to 
	*/
	if( remapBL ) 
	  {
	    neighborZ = nbz; 
	    for(int i = 0; i < chain->traln->getNumBranches(); ++i)
	      traln->setBranchLengthSave(traln->getBranchLength(nnb, i), i,nb); 
	    
	    assert(chain->traln->getNumBranches() == 1 ); 
	    swpDouble(rec->neighborBls,rec->nextNeighborBls ); 
	    rec->pruningBranch = invertBranch(rec->pruningBranch); 
	  }
	
	
	insertNodeIntoBranch( chain, 
			      constructBranch(prunedSubtreePtr->number, prunedSubtreePtr->back->number ),
			      rec->insertBranch, 
			      insertBranchPtr->z,
			      neighborZ); 	
      }
      break; 
    default: assert(0); 
    }

  debug_checkTreeConsistency(chain->traln->getTr()); 

  assert(branchExists(tr, constructBranch(rec->pruned.thatNode, rec->insertBranch.thisNode ))); 
  assert(branchExists(tr, constructBranch(rec->pruned.thatNode, rec->insertBranch.thatNode ))); 
} 



static void extended_spr_reset(Chain * chain, proposalFunction *pf)
{
  tree *tr = chain->traln->getTr(); 
  topoRecord
    *topoRec = pf->remembrance.topoRec; 
  

  /* HACK to find the other node in the subtree */  
  nodeptr p = tr->nodep[topoRec->pruned.thatNode] ; 
  while(p->back->number == topoRec->insertBranch.thisNode 
	|| p->back->number == topoRec->insertBranch.thatNode )
    p = p->next; 
  p = p->back; 
  assert(p->number != topoRec->insertBranch.thisNode  && p->number != topoRec->insertBranch.thatNode); 
  /* END */

  pruneBranch(chain, constructBranch(topoRec->pruned.thatNode, p->number) , topoRec->bls);

  insertNodeIntoBranch(chain, 
		       constructBranch(topoRec->pruned.thatNode, p->number), 
		       topoRec->pruningBranch, 
		       topoRec->neighborBls, topoRec->nextNeighborBls) ; 

#ifdef DEBUG_SHOW_TOPO_CHANGES
  debug_printNodeEnvironment(chain,  topoRec->pruned.thatNode ); 
  debug_printNodeEnvironment(chain,   topoRec->pruningBranch.thisNode ); 
  debug_printNodeEnvironment(chain,   topoRec->pruningBranch.thatNode ); 
  debug_printNodeEnvironment(chain,   topoRec->insertBranch.thisNode ); 
  debug_printNodeEnvironment(chain,   topoRec->insertBranch.thatNode ); 
#endif

  debug_checkTreeConsistency(chain->traln->getTr() );
  debug_printTree(chain);
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

  int model = drawRandInt(chain, chain->traln->getNumberOfPartitions());  
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





//NOTE: should only be called at the very beginning. Afterwards a probability sum of 1.0 is not required.


#if  0
void normalizeProposalWeights(Chain *chain)
{
  double sum = 0 ; 
  for(int i = 0; i < NUM_PROPOSALS;++i)
    sum += chain->proposalWeights[i]; 
  
  for(int i = 0; i < NUM_PROPOSALS;++i)
    chain->proposalWeights[i] /= sum; 
  
  findLogisticT(chain);
}
#endif



static void simple_model_proposal_apply(Chain *chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 

  int model = drawRandInt(chain, chain->traln->getNumberOfPartitions());
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
	      curv = traln->getSubstRate(model, changeState); 
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

  int model = drawRandInt(chain, chain->traln->getNumberOfPartitions());
  perPartitionInfo *info = pf->remembrance.partInfo;

  pInfo *partition = chain->traln->getPartition( model); 
  int numRates = (partition->states * partition->states - partition->states) / 2; 

  recordSubsRates(chain, model , numRates, info->substRates);

  info->modelNum = model;
  info->numRates = numRates; 

  double r[numRates];

  int* list  = (int*)exa_calloc(numRates, sizeof(int)); 

     
  drawDirichletExpected(chain, r, partition->substRates, 50.0, numRates);
  chain->hastings=densityDirichlet(partition->substRates, r, numRates) / densityDirichlet(r, partition->substRates, numRates);    
  
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

  int model = drawRandInt(chain,chain->traln->getNumberOfPartitions());; 

  pInfo *partition = chain->traln->getPartition( model) ; 

  perPartitionInfo *info = pf->remembrance.partInfo; 

  int numRates  = (partition->states  * partition->states -partition->states ) / 2 ;  
  info->numRates = numRates; 

  recordSubsRates(chain, model, numRates, info->substRates);
  int state, randNumber;
  double new_value,curv;
  double r;
  
  randNumber=drawRandInt(chain,numRates);
  int *perm = (int*)exa_calloc(numRates, sizeof(int)); 
  // int perm[numRates];
  drawPermutation(chain,perm, numRates);

  for(state = 0;state < randNumber ; state ++)
    {           
      curv = traln->getSubstRate(model, perm[state]);;
      r =  drawRandBiUnif(chain,curv);

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
  int model = drawRandInt(chain,chain->traln->getNumberOfPartitions()); 

  perPartitionInfo *info =  pf->remembrance.partInfo; 
  
  pInfo *partition = chain->traln->getPartition(model) ; 

  info->modelNum = model; 
  info->numRates = partition->states ;   
  int numRates = info->numRates; 
  recordSubsRates(chain, model, info->numRates, info->substRates);
  //choose a random set parameter,
  //with uniform probabilities

  int randState= drawRandInt(chain,numRates);

  double new_value,curv;
  double r;
  
  //int state=drawRandInt(chain->modelRemem.numSubsRates);
  

  curv = traln->getSubstRate(model, randState);
  r =  drawRandBiUnif(chain,curv);

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
  
  int model = drawRandInt(chain,chain->traln->getNumberOfPartitions());

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
      r =  drawRandBiUnif(chain,curv);

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

  assert(partition->dataType = DNA_DATA);
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
      r = drawRandDouble01(chain);

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
      r = drawRandBiUnif(chain, real_z);//draw from [real_z/2,2*real_z]

    
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
       r = drawRandExp(chain, 1.0/real_z);

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
  int target_branch = drawRandInt(chain,(chain->traln->getTr()->mxtips * 2) - 3); 

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
  pushStack(pf->remembrance.modifiedPath, b); 
  pf->remembrance.modifiedPath->content[0].length[0] = z; 
  
  /* pf->remembrance.topoRec->whichBranch = target_branch;  */
  /* pf->remembrance.insertBranch = constructBranch(p->) */

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

static void random_branch_length_proposal_reset(Chain * chain, proposalFunction *pf)
{
  /* branch b = pf->remembrance.modifiedPath->content[0];  */
  branch b = popStack(pf->remembrance.modifiedPath); 
  nodeptr p = findNodeFromBranch(chain->traln->getTr(), b); 
  chain->traln->setBranchLengthSave(b.length[0],0,p); 

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

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    chain->traln->setFrequencySave(prevFrequRates[i], model, i); 

  chain->traln->initRevMat(model);
}

static void recordFrequRates(Chain *chain, int model, int numFrequRates, double *prevFrequRates)
{
  pInfo *partition = chain->traln->getPartition(model);

  assert(partition->dataType = DNA_DATA);
  int i;
  for(i=0; i<numFrequRates; i++)
    prevFrequRates[i] = chain->traln->getFrequency(model,i);
}



static void frequencySliderApply(Chain *chain, proposalFunction *pf)
{  
  int model = drawRandInt(chain, chain->traln->getNumberOfPartitions()); 
  perPartitionInfo *info = pf->remembrance.partInfo; 
  pInfo*partition = chain->traln->getPartition(model); 

  int numFreq = partition->states; 
  info->modelNum = model; 
  info->numFreqs = numFreq; 

  recordFrequRates(chain, model, numFreq, info->frequencies); 
  
  int paramToChange = drawRandInt(chain, numFreq); 
  double curv = chain->traln->getFrequency(model, paramToChange);   
  double newVal = fabs(drawFromSlidingWindow(chain,curv, pf->parameters.slidWinSize)); 
  chain->traln->setFrequencySave(newVal,model, paramToChange); 

  double sum = 0;
  for(int i = 0; i < numFreq; ++i)
    sum += chain->traln->getFrequency(model,i); 
  for(int i = 0; i < numFreq; ++i)
    chain->traln->setFrequencySave(
				   chain->traln->getFrequency(model,i) / sum,
				   model,i); 

  
  chain->traln->initRevMat(model);
  // exa_initReversibleGTR(chain, model);
}



void frequency_proposal_apply(Chain * chain, proposalFunction *pf)
{
  int model  = drawRandInt(chain,chain->traln->getNumberOfPartitions());
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
      r[state] = drawRandBiUnif(chain,curv);       
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
   
  int model  = drawRandInt(chain,chain->traln->getNumberOfPartitions());
  perPartitionInfo *info = pf->remembrance.partInfo; 

  pInfo *partition = chain->traln->getPartition( model); 

  int numFreq = partition->states; 
  info->modelNum = model;   
  info->numFreqs = numFreq; 

  recordFrequRates(chain, model, numFreq, info->frequencies);
  
  double* r  =  (double*)exa_calloc(numFreq, sizeof(double)); 
  
  
        
  //drawRandDirichlet(chain, r, partition->frequencies, 1.0, numFreq);
  drawDirichletExpected(chain, r, partition->frequencies, 50.0, numFreq);
  chain->hastings=densityDirichlet(partition->frequencies, r, numFreq) / densityDirichlet(r, partition->frequencies, numFreq);

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

  assert(partition->dataType = DNA_DATA);
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
  branch b =  drawBranchUniform(chain); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  clearStack(pf->remembrance.modifiedPath); 
  pushStack(pf->remembrance.modifiedPath, b);
  pf->remembrance.modifiedPath->content[0].length[0]  = chain->traln->getBranchLength( p,0); 

  double zOld = chain->traln->getBranchLength( p,0); 
  double win =  pf->parameters.slidWinSize ; 
  double realZ = branchLengthToReal(tr, zOld); 
  double blNew = branchLengthToInternal(tr,fabs(drawFromSlidingWindow(chain, realZ, win))); 
  chain->traln->setBranchLengthSave(blNew, 0,p);
  // p->z[0] = p->back->z[0] = blNew; 
}

static void branchLengthMultiplierApply(Chain *chain, proposalFunction *pf)
{
  auto traln =  chain->traln; 
  tree *tr = chain->traln->getTr(); 
  branch b =  drawBranchUniform(chain); 

  nodeptr p = findNodeFromBranch(tr, b); 

  clearStack(pf->remembrance.modifiedPath); 
  pushStack(pf->remembrance.modifiedPath, b);
  assert(stackLength(pf->remembrance.modifiedPath) == 1 ); 
  pf->remembrance.modifiedPath->content[0].length[0]  = traln->getBranchLength( p,0); 

  double
    multiplier = drawMultiplier(chain, pf->parameters.multiplier); 
  assert(multiplier > 0.); 
  
  /* TODO how do we do that wiht multiple bls per branch?  */
  assert(chain->traln->getNumBranches() == 1); 

  double newZ = pow( traln->getBranchLength( p,0), multiplier) ; 
  /* printf("\t\tBL_MULTI: %g * %g => %g\n", branchLengthToReal(tr, p->z[0]), multiplier, branchLengthToReal(tr, newZ));  */

  /* according to lakner2008  */
  chain->hastings *= multiplier; 
 
  /* just doing it for one right here */
  traln->setBranchLengthSave(newZ, 0, p); 
}


/**
   @brief tunes the BL multiplier
 */ 
static void autotuneMultiplier(Chain *chain, proposalFunction *pf)
{
  double *parameter = &(pf->parameters.multiplier); 

  // successCtr *ctr = &(pf->sCtr); 
  SuccessCtr *ctr = &(pf->sCtr); 

  int batch = chain->currentGeneration  / gAInfo.tuneFreq; 

  double newParam = tuneParameter(batch, ctr->getRatioInLastInterval(), *parameter, FALSE); 

#ifdef DEBUG_PRINT_TUNE_INFO
  printInfo(chain, "%s\tratio=%f\t => %s %f to %f\n", pf->name, getRatioLocal(ctr), (newParam < *parameter ) ? "reducing" : "increasing", *parameter, newParam);
#endif

  *parameter = newParam; 
  ctr->reset();   
}


/**
   @brief autotunes sliding windows 
   
*/ 
static void autotuneSlidingWindow(Chain *chain, proposalFunction *pf)
{
  double *parameter = &(pf->parameters.slidWinSize); 
  SuccessCtr *ctr = &(pf->sCtr); 
  double newParam = tuneParameter(chain->currentGeneration / gAInfo.tuneFreq,
				  ctr->getRatioInLastInterval(), 
				  *parameter, FALSE  ); 
  
#ifdef DEBUG_PRINT_TUNE_INFO
  printInfo(chain, "%s\tratio=%f\t => %s %f to %f\n", pf->name, getRatioLocal(ctr), (newParam < *parameter ) ? "reducing" : "increasing", *parameter, newParam);
#endif

  *parameter = newParam; 
  ctr->reset(); 
}



// TODO this below did not work at all 
/**
   @brief autotunes the stop probability of extended topological moves 
*/ 
// static void autotuneStopProp(Chain *chain, proposalFunction *pf) 
// {
//   const double minimum = 0.01; 
//   const double maximum = 0.95; 

//   double *parameter = &(pf->parameters.eSprStopProb); 
//   successCtr *ctr = &(pf->sCtr); 
//   double newParam = tuneParameter( chain->currentGeneration / gAInfo.tuneFreq, 
// 				   getRatioLocal(&(pf->sCtr)),
// 				   *parameter, 
// 				   TRUE ); 

// #ifdef DEBUG_PRINT_TUNE_INFO
//   printInfo(chain, "%s\tratio=%f\t => %s %f to %f\n", pf->name, getRatioLocal(ctr), (newParam < *parameter ) ? "reducing" : "increasing", *parameter, newParam);
// #endif
  
//   *parameter = fmax(minimum,fmin(newParam,maximum)); 
//   resetCtr(ctr);     
// }


void branchLengthReset(Chain *chain, proposalFunction *pf)
{
  branch b = popStack(pf->remembrance.modifiedPath); 
  tree *tr = chain->traln->getTr(); 
  nodeptr p = findNodeFromBranch(tr, b); 
  chain->traln->setBranchLengthSave(b.length[0], 0,p); 
  // p->z[0] = p->back->z[0] = b.length[0]; 
}




void resetGammaMulti(Chain *chain, proposalFunction *pf)
{
  int model = pf->remembrance.partInfo->modelNum; 
  chain->traln->setAlphaSave(pf->remembrance.partInfo->alpha, model); 
  chain->traln->discretizeGamma(model); 
}

 
void applyGammaMultiplier(Chain *chain, proposalFunction *pf)
{
  double multi = drawMultiplier(chain, pf->parameters.multiplier);   
  int model = drawRandInt(chain, chain->traln->getNumberOfPartitions());
  pInfo *partition = chain->traln->getPartition( model); 

  
  perPartitionInfo *pinfo = pf->remembrance.partInfo; 
  pinfo->modelNum = model;   
  pinfo->alpha = chain->traln->getAlpha(model);
  
  chain->traln->setAlphaSave(partition->alpha * multi, model); 
  chain->traln->discretizeGamma(model);
}


 
void initProposalFunction( proposal_type type, initParamStruct *initParams, proposalFunction **result)
{
  if(initParams->initWeights[type] == 0)
    {
      *result = NULL; 
      return ; 
    }  

  *result = (proposalFunction*)exa_calloc(1,sizeof(proposalFunction));   

  proposalFunction *ptr = *result; 
  ptr->ptype = (proposal_type)type; 
  ptr->initWeight = initParams->initWeights[type]; 
  ptr->currentWeight = ptr->initWeight; 


  /* 
     TODO@kassian
     some proposals may be broken      
  */

  switch(type)
    {
    case GUIDED_SPR:
      ptr->name = "guidedSpr"; 
      ptr->apply_func = applyGuidedSPR; 
      ptr->eval_lnl = evalGuidedSPR; 
      ptr->reset_func = resetGuidedSPR; 
      ptr->remembrance.modifiedPath = NULL; 
      createStack(&(ptr->remembrance.modifiedPath)); 
      ptr->category = TOPOLOGY; 
      ptr->parameters.radius =  initParams->initGuidedSPR; 
      ptr->autotune = NULL;       
      break; 
    case E_SPR:
      ptr->eval_lnl = sprEval; 
      /* ptr->autotune = autotuneStopProp; */
      ptr->remembrance.modifiedPath = (path*)exa_calloc(1,sizeof(path)); 
      ptr->apply_func = applyExtendedSPR; 
      ptr->reset_func = resetESPR; 
      ptr->name = "eSPR"; 
      ptr->param2.multiplier = INIT_ESPR_MULT; 
      ptr->parameters.eSprStopProb = initParams->eSprStopProb; 
      ptr->category = TOPOLOGY; 
      break; 
    case UPDATE_MODEL: 	
      ptr->eval_lnl = onePartitionEval; 
      ptr->autotune = autotuneSlidingWindow; 
      ptr->apply_func = simple_model_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->parameters.slidWinSize = INIT_RATE_SLID_WIN;
      ptr->category = SUBSTITUTION_RATES; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "modelSlidWin"; 
      break; 
    case UPDATE_GAMMA:      	
      ptr->eval_lnl = onePartitionEval;
      ptr->autotune = autotuneSlidingWindow; 
      ptr->apply_func = simple_gamma_proposal_apply; 
      ptr->reset_func = simple_gamma_proposal_reset; 
      ptr->parameters.slidWinSize = INIT_RATE_SLID_WIN; 
      ptr->category = RATE_HETEROGENEITY; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "gammaSlidWin"; 
      break; 
    case UPDATE_GAMMA_EXP: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = simple_gamma_proposal_apply; 
      ptr->reset_func = simple_gamma_proposal_reset; 
      ptr->category = RATE_HETEROGENEITY; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "gammaExp"; 
      break; 
    case UPDATE_SINGLE_BL: 
      ptr->eval_lnl = evalBranch;
      ptr->autotune = autotuneSlidingWindow; 
      ptr->remembrance.modifiedPath = NULL; 
      createStack(&(ptr->remembrance.modifiedPath)); 
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
    case UPDATE_MODEL_DIRICHLET:
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = model_dirichlet_proposal_apply; 
      ptr->reset_func = simple_model_proposal_reset; 
      ptr->name = "modelDirichlet"; 
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
    case UPDATE_FREQUENCIES_DIRICHLET: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = frequency_dirichlet_proposal_apply; 
      ptr->reset_func = frequency_proposal_reset; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "freqDirichlet"; 
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
      ptr->remembrance.modifiedPath =  NULL; 
      createStack(&(ptr->remembrance.modifiedPath));
      ptr->parameters.multiplier = INIT_BL_MULT; 
      ptr->name = "branchMult"; 
      ptr->category = BRANCH_LENGTHS; 
      break; 
    case FREQUENCY_SLIDER: 
      ptr->eval_lnl = onePartitionEval; 
      ptr->apply_func = frequencySliderApply; 
      ptr->autotune = autotuneSlidingWindow; 
      ptr->reset_func = frequency_proposal_reset; 
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name = "freqSlider"; 
      ptr->category = FREQUENCIES; 
      ptr->parameters.slidWinSize = INIT_FREQ_SLID_WIN; 
      break; 
    case ST_NNI: 
      ptr->eval_lnl = eval_st_nni; 
      ptr->apply_func = apply_st_nni; 
      ptr->autotune = NULL; 
      ptr->reset_func = reset_st_nni; 
      ptr->remembrance.modifiedPath = NULL; 
      createStack(&(ptr->remembrance.modifiedPath)); 
      ptr->name = "stNNI"; 
      ptr->category = TOPOLOGY; 
      ptr->parameters.multiplier = INIT_NNI_MULT; 
      break;       
    case TL_MULT: 
      ptr->eval_lnl = dummy_eval;
      ptr->apply_func = applyTLMult; 
      ptr->reset_func = resetTLMult; 
      ptr->autotune = autotuneMultiplier; 
      ptr->remembrance.multiplier = 0; 
      ptr->name = "tl-mult"; 
      ptr->category = BRANCH_LENGTHS; 
      ptr->parameters.multiplier = INIT_TL_MULTI; 
      break;       
    case NODE_SLIDER: 
      ptr->apply_func = applyNodeSlider; 
      ptr->eval_lnl = evaluateNodeSlider; 
      ptr->reset_func = resetNodeSlider;  
      ptr->autotune = autotuneMultiplier; 
      ptr->remembrance.modifiedPath = NULL; 
      createStack(&(ptr->remembrance.modifiedPath)); 
      ptr->name = "nodeSlider"; 
      ptr->category = BRANCH_LENGTHS; 
      ptr->parameters.multiplier = INIT_NODE_SLIDER_MULT; 
      break; 
    case GAMMA_MULTI: 
      ptr->apply_func = applyGammaMultiplier;       
      ptr->eval_lnl = onePartitionEval; 
      ptr->reset_func = resetGammaMulti; 
      ptr->autotune = autotuneMultiplier;  
      ptr->remembrance.partInfo = (perPartitionInfo*)exa_calloc(1,sizeof(perPartitionInfo)); 
      ptr->name =  "gammaMulti"; 
      ptr->category = RATE_HETEROGENEITY; 
      ptr->parameters.multiplier = INIT_GAMMA_MULTI;       
      break; 
      /* TODO re-install PROPOSALADD anchor for script   */
    default : 
      {
	printf("unknown value %d\n",type ); 
	assert(0) ; 
      }
    }  
}

