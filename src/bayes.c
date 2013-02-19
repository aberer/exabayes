#include "common.h"

#include "axml.h"
#include "proposalStructs.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"

#include "proposals.h"



/* TODO commented this out, since there are some problems with it,
   when we build with the PLL */
/* #define WITH_PERFORMANCE_MEASUREMENTS */


#define _USE_NCL_PARSER


void makeRandomTree(tree *tr); 

#ifdef _USE_NCL_PARSER
#include "nclConfigReader.h"
void addInitParameters(state *curstate, initParamStruct *initParams)
{
  curstate->proposalWeights[E_SPR] = initParams->initSPRWeight; 
  curstate->proposalWeights[UPDATE_GAMMA] = initParams->initGammaWeight;
  curstate->proposalWeights[UPDATE_MODEL] = initParams->initModelWeight; 
  curstate->proposalWeights[UPDATE_SINGLE_BL] = initParams->initSingleBranchWeight; 
  curstate->proposalWeights[UPDATE_SINGLE_BL_EXP] = initParams->initSingleBranchExpWeight; 
  curstate->numGen = initParams->numGen; 
  curstate->penaltyFactor = initParams->initPenaltyFactor; 
  curstate->samplingFrequency = initParams->samplingFrequency; 
  curstate->eSprStopProb = initParams->eSprStopProb; 
}
#else 
int parseConfig(state *theState);
#endif



static void traverse_branches_set_fixed(nodeptr p, int *count, state * s, double z )
{
  nodeptr q;
  int i;
  
  for( i = 0; i < s->tr->numBranches; i++)
    p->z[i] = p->back->z[i] = z;  
  *count += 1;


  if (! isTip(p->number, s->tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverse_branches_set_fixed(q->back, count, s, z);
	  q = q->next;
	}   
      newviewGeneric(s->tr, p, FALSE);     // not sure if we need this
    }
}


void initDefaultValues(state *theState, tree *tr)
{
  theState->curprior = 1; 
  theState->hastings = 1; 
  theState->currentGeneration = 0; 

  theState->brLenRemem.bl_sliding_window_w = 0.005;
  theState->brLenRemem.bl_prior = 1.0;
  theState->brLenRemem.bl_prior_exp_lambda = 0.1 ;
  //this can be extended to more than one partition, but one for now
  theState->modelRemem.model = 0;

  theState->modelRemem.rt_sliding_window_w = 0.5;
  theState->modelRemem.nstates = tr->partitionData[theState->modelRemem.model].states; /* 4 for DNA */
  theState->modelRemem.numSubsRates = (theState->modelRemem.nstates * theState->modelRemem.nstates - theState->modelRemem.nstates) / 2; /* 6 for DNA */
  theState->modelRemem.curSubsRates = (double *) malloc(theState->modelRemem.numSubsRates * sizeof(double));
  theState->gammaRemem.gm_sliding_window_w = 0.75;
  theState->brLenRemem.single_bl_branch = -1;


  theState->proposalWeights[E_SPR] = 0.0; 
  theState->proposalWeights[UPDATE_MODEL] = 0.0; 
  theState->proposalWeights[UPDATE_GAMMA] = 0.0; 
  theState->proposalWeights[UPDATE_SINGLE_BL] = 0.0;   
  theState->proposalWeights[UPDATE_SINGLE_BL_EXP] = 0.0;   
  
  theState->numGen = 1000000;
  theState->penaltyFactor = 0.0;
}


void readConfig(state *curstate, tree *tr)
{
  initDefaultValues(curstate, tr);
#ifdef _USE_NCL_PARSER
  initParamStruct *initParams = NULL; 
  parseConfigWithNcl(configFileName, &initParams);   
  addInitParameters(curstate, initParams); 
#else  
  parseConfig(curstate); 
#endif
  normalizeProposalWeights(curstate); 
  
}


state *state_init(tree *tr, analdef * adef)
{
  state *curstate  =(state *)calloc(1,sizeof(state));

  nodeptr *list = (nodeptr *)malloc(sizeof(nodeptr) * 2 * tr->mxtips);
  curstate->list = list;

  curstate->tr = tr;

  curstate->modelRemem.adef = adef;

  assert(curstate != NULL);

  return curstate;
}


static node *find_tip( node *n, tree *tr ) {
  if( isTip(n->number, tr->mxtips) ) {
    return n;
  } else {
    return find_tip( n->back, tr );
  }
  
}

#if 0
static char *Tree2StringRecomREC(char *treestr, tree *tr, nodeptr q, boolean printBranchLengths)
{
  char  *nameptr;            
  double z;
  nodeptr p = q;

  if(isTip(p->number, tr->mxtips)) 
    {               
      nameptr = tr->nameList[p->number];     
      sprintf(treestr, "%s", nameptr);
      while (*treestr) treestr++;
    }
  else 
    {                      
      while(!p->x)
	p = p->next;
      *treestr++ = '(';
      treestr = Tree2StringRecomREC(treestr, tr, q->next->back, printBranchLengths);
      *treestr++ = ',';
      treestr = Tree2StringRecomREC(treestr, tr, q->next->next->back, printBranchLengths);
      if(q == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2StringRecomREC(treestr, tr, q->back, printBranchLengths);
	}
      *treestr++ = ')';                    
      // write innernode as nodenum_b_nodenumback
#if 0
      sprintf(treestr, "%d", q->number);
      while (*treestr) treestr++;
      *treestr++ = 'b';                    
      sprintf(treestr, "%d", p->back->number);
      while (*treestr) treestr++;
#endif
    
    }

  if(q == tr->start->back) 
    {              
      if(printBranchLengths)
	sprintf(treestr, ":0.0;\n");
      else
	sprintf(treestr, ";\n");                  
    }
  else 
    {                   
      if(printBranchLengths)          
	{
	  //sprintf(treestr, ":%8.20f", getBranchLength(tr, SUMMARIZE_LH, p));                 
	  assert(tr->fracchange != -1.0);
	  z = q->z[0];
	  if (z < zmin) 
	    z = zmin;        
	  sprintf(treestr, ":%8.20f", -log(z) * tr->fracchange);               
	}
      else            
	sprintf(treestr, "%s", "\0");         
    }

  while (*treestr) treestr++;
  return  treestr;
}
#endif




void mcmc(tree *tr, analdef *adef)
{    
  initRNG(seed);
  
  tr->start = find_tip(tr->start, tr );

  assert( isTip(tr->start->number, tr->mxtips ));
  
  state *curstate = state_init(tr, adef); 

  readConfig(curstate, tr);
  
  if(processID == 0 )
    initializeOutputFiles(curstate);

  int count = 0;
  traverse_branches_set_fixed( tr->start, &count, curstate, 0.65 );

  /* makeRandomTree(tr); */

  evaluateGeneric(tr, tr->start, TRUE);
  PRINT( "after reset start: %f\n\n", tr->likelihood );


  /* beginning of the MCMC chain */
  while(curstate->currentGeneration < curstate->numGen)
    step(curstate);

  if(processID == 0)
    finalizeOutputFiles(curstate);
}


