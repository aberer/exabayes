#include "common.h"

#include "axml.h"
#include "proposalStructs.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"

#include "proposals.h"


#define _USE_NCL_PARSER


void makeRandomTree(tree *tr); 



#ifdef _USE_NCL_PARSER
#include "nclConfigReader.h"
void addInitParameters(state *curstate, initParamStruct *initParams)
{
  curstate->proposalWeights[E_SPR] = initParams->initSPRWeight; 
  curstate->proposalWeights[UPDATE_GAMMA] = initParams->initGammaWeight;
  curstate->proposalWeights[UPDATE_GAMMA_EXP] = initParams->initGammaExpWeight;
  curstate->proposalWeights[UPDATE_MODEL] = initParams->initModelWeight; 
  curstate->proposalWeights[UPDATE_SINGLE_BL] = initParams->initSingleBranchWeight; 
  curstate->proposalWeights[UPDATE_SINGLE_BL_EXP] = initParams->initSingleBranchExpWeight; 
curstate->proposalWeights[UPDATE_SINGLE_BL_BIUNIF] = initParams->initSingleBranchBiunifWeight;
  //PROPOSALADD addInitParameters NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
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
  //printf("current BL at %db%d: %f\n", p->number, p->back->number, p->z[0]);
  
  for( i = 0; i < s->tr->numBranches; i++)
    {
   
      p->z[i] = p->back->z[i] = z;
    }
  
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


void initDefaultValues(state *theState)
{
  
  theState->proposalWeights[E_SPR] = 0.0; 
  theState->proposalWeights[UPDATE_MODEL] = 0.0; 
  theState->proposalWeights[UPDATE_GAMMA] = 0.0; 
  theState->proposalWeights[UPDATE_GAMMA_EXP] = 0.0; 
  theState->proposalWeights[UPDATE_SINGLE_BL] = 0.0;   
  theState->proposalWeights[UPDATE_SINGLE_BL_EXP] = 0.0;   
theState->proposalWeights[UPDATE_SINGLE_BL_BIUNIF] = 0.0;
  //PROPOSALADD initDefaultValues NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
  
  theState->numGen = 1000000;
  theState->penaltyFactor = 0.0;
}




/* TODO commented this out, since there are some problems with it,
   when we build with the PLL */
/* #define WITH_PERFORMANCE_MEASUREMENTS */





void readConfig(state *curstate)
{
  initDefaultValues(curstate);
#ifdef _USE_NCL_PARSER
  initParamStruct *initParams = NULL; 
  parseConfigWithNcl(configFileName, &initParams);   
  addInitParameters(curstate, initParams); 
#else  
  parseConfig(curstate); 
#endif
  normalizeProposalWeights(curstate); 
  
}


state *state_init(tree *tr, analdef * adef, double bl_w, double rt_w, double gm_w, double bl_p)
{
  state *curstate  =(state *)calloc(1,sizeof(state));
  nodeptr *list = (nodeptr *)malloc(sizeof(nodeptr) * 2 * tr->mxtips);
  curstate->list = list;

  curstate->tr = tr;
  curstate->brLenRemem.bl_sliding_window_w = bl_w;
  curstate->brLenRemem.bl_prior = 1.0;
  curstate->brLenRemem.bl_prior_exp_lambda = bl_p;
  //this can be extended to more than one partition, but one for now
  curstate->modelRemem.model = 0;
  curstate->modelRemem.adef = adef;
  curstate->modelRemem.rt_sliding_window_w = rt_w;
  curstate->modelRemem.nstates = tr->partitionData[curstate->modelRemem.model].states; /* 4 for DNA */
  curstate->modelRemem.numSubsRates = (curstate->modelRemem.nstates * curstate->modelRemem.nstates - curstate->modelRemem.nstates) / 2; /* 6 for DNA */
  curstate->modelRemem.curSubsRates = (double *) malloc(curstate->modelRemem.numSubsRates * sizeof(double));
  curstate->gammaRemem.gm_sliding_window_w = gm_w;
  assert(curstate != NULL);
  
  curstate->brLenRemem.single_bl_branch = -1;

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
  //allocate states
  double bl_prior_exp_lambda = 0.1;
  double bl_sliding_window_w = 0.005;
  double gm_sliding_window_w = 0.75;
  double rt_sliding_window_w = 0.5;

  initRNG(seed);
  
  tr->start = find_tip(tr->start, tr );

  assert( isTip(tr->start->number, tr->mxtips ));
  
  state *curstate = state_init(tr, adef,  bl_sliding_window_w, rt_sliding_window_w, gm_sliding_window_w, bl_prior_exp_lambda); 

  readConfig(curstate);

  for(int prop=0; prop<NUM_PROPOSALS;prop++)
  {
  printf("%f ",curstate->proposalWeights[prop]);
  }
  printf("\n");
  
  if(processID == 0 )
    initializeOutputFiles(curstate);

  int count = 0;
  traverse_branches_set_fixed( tr->start, &count, curstate, 0.65 );

  makeRandomTree(tr);

  evaluateGeneric(tr, tr->start, TRUE);
  PRINT( "after reset start: %f\n\n", tr->likelihood );
  
  curstate->curprior = 1;
  curstate->hastings = 1;
  curstate->currentGeneration = 0; 

#ifdef WITH_PERFORMANCE_MEASUREMENTS  
  perf_timer all_timer = perf_timer_make();
#endif
  
  /* beginning of the MCMC chain */
  while(curstate->currentGeneration < curstate->numGen)
    step(curstate);

  if(processID == 0)
    finalizeOutputFiles(curstate);
}


