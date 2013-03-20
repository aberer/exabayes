/**
   @file output.c
   
   @brief Functions printing to console or files. 
 */ 


#include "common.h"
#include "config.h"
#include "axml.h"
#include "proposalStructs.h"
#include "main-common.h"
#include "globals.h"


#include "chain.h"
#include "output.h"

#include "adapters.h"		/*  */

extern double masterTime; 

void makeFileNames()
{
  int infoFileExists = 0;

  strcpy(infoFileName, workdir); 
  strcat(infoFileName, PROGRAM_NAME);
  strcat(infoFileName, "_infoFile."); 
  strcat(infoFileName, run_id); 

  infoFileExists = filexists(infoFileName);

  if(infoFileExists)
  {
    printf("%s output files with the run ID <%s> already exist \n", PROGRAM_NAME, run_id);
    printf("in directory %s ...... exiting\n", workdir);
#ifdef PRODUCTIVE
    exit(-1);
#endif
  }
}



double exabayes_getBranchLength(state *chain, int perGene, nodeptr p)
{
  tree *tr = chain->tr; 

  double 
    z = 0.0,
    x = 0.0;

  assert(perGene != NO_BRANCHES);
	      
  if(NOT hasPergeneBL(tr))
    {
      assert(tr->fracchange != -1.0);
      z = p->z[0];
      if (z < zmin) 
	z = zmin;      	 
      
      x = -log(z) * tr->fracchange; 
    }
  else
    {
      if(perGene == SUMMARIZE_LH)
	{
	  int 
	    i;
	  
	  double 
	    avgX = 0.0;
		      
	  for(i = 0; i < getNumBranches(tr); i++)
	    {

	      assert(getPcontr(chain,i) != -1.0);
	      assert(getFracChange(chain,i) != -1.0);

	      z = p->z[i];
	      if(z < zmin) 
		z = zmin;      	 
	      x = -log(z) * getFracChange(chain,i);
	      avgX += x * getPcontr(chain,i);
	    }

	  x = avgX;
	}
      else
	{	
	  assert(getFracChange(chain,perGene) != -1.0);
	  assert(perGene >= 0 && perGene < getNumBranches(tr));
	  
	  z = p->z[perGene];
	  
	  if(z < zmin) 
	    z = zmin;      	 
	  
	  x = -log(z) * getFracChange(chain,perGene);
	}
    }

  return x;
}


char *Tree2stringNexus(char *treestr, state *chain , nodeptr p, int perGene )
{      
  tree *tr = chain->tr; 

  if(isTip(p->number, tr->mxtips)) 
    {	       	        
      sprintf(treestr, "%d", p->number); 
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = Tree2stringNexus(treestr, chain, p->next->back, perGene);
      *treestr++ = ',';
      treestr = Tree2stringNexus(treestr, chain, p->next->next->back, perGene);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2stringNexus(treestr, chain, p->back, perGene);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back) 
    {
      sprintf(treestr, ";"); 
      return treestr;
    }

  sprintf(treestr, ":%8.20f", exabayes_getBranchLength(chain, perGene, p));
  
  while (*treestr) treestr++;
  return  treestr;
}


/* TODO maybe generate random run id */


static void printNexusTreeFileStart(state *curstate)
{
  FILE *fh = curstate->topologyFile; 

  fprintf(fh, "#NEXUS\n[ID: %s]\n[Param: tree]\nbegin trees;\n\ttranslate\n", "TODO" ); 

  for(int i = 0; i < curstate->tr->mxtips-1; ++i)
    fprintf(fh, "\t\t%d %s,\n", i+1, curstate->tr->nameList[i+1]); 
  fprintf(fh, "\t\t%d %s;\n", curstate->tr->mxtips, curstate->tr->nameList[curstate->tr->mxtips]);
  /* fclose(fh);  */
}


static void printParamFileStart(state *curstate)
{
  /* FILE *fh = myfopen(outputParamFile, "w");  */
  FILE *fh = curstate->outputParamFile; 
  
  char *tmp = "TODO"; 
  fprintf(fh, "[ID: %s]\n", tmp);
  fprintf(fh, "Gen\tLnL\tTL\
\tr(A<->C)\tr(A<->G)\tr(A<->T)\tr(C<->G)\tr(C<->T)\tr(G<->T)\t\
pi(A)\tpi(C)\tpi(G)\tpi(T)\
\talpha\tpinvar\n"); 
  fflush(fh); 
}

static void printParams(state *curstate)
{
  int model = 0; 

  pInfo *partition = getPartition(curstate, model); 

  /* FILE *fh = myfopen(outputParamFile, "a");  */
  FILE *fh = curstate->outputParamFile; 

  /* TODO tree length */
  fprintf(fh, "%d\t%f\t%f", curstate->currentGeneration,curstate->tr->likelihood, -1. ); 

  /* TODO what about multiple models?  */
  /* TODO that will not work in the future ...  */

  for(int i = 0; i< curstate->modelRemem.numSubsRates;  ++i )
    fprintf(fh, "\t%f" , partition->substRates[i]); 

  
  /* TODO works on dna only currently  */
  for(int i = 0; i < 4; ++i)
    fprintf(fh, "\t%f", partition->frequencies[i]); 

  /* TODO what is pinvar?  */
  fprintf(fh, "\t%f\t%f", partition->alpha, -1.  ); 

  fprintf(fh, "\n");

  fflush(fh); 
}



void initializeOutputFiles(state *curstate)
{
  printNexusTreeFileStart(curstate); 
  printParamFileStart(curstate);
}


void finalizeOutputFiles(state *curstate)
{
  FILE *fh = curstate->topologyFile; 

  /* topo file  */
  fprintf(fh, "end;\n"); 
  fclose(fh); 
}



  /* TODO what about per model brach lengths? how does mrB do this? */
static void exabayes_printTopology(state *curstate)
{
  /* FILE *fh = myfopen(topologyFile, "a");   */
  FILE *fh = curstate->topologyFile; 
  memset(curstate->tr->tree_string, 0, curstate->tr->treeStringLength * sizeof(char) ); 
  
  Tree2stringNexus(curstate->tr->tree_string, curstate,  curstate->tr->start->back, 0 ); 
  fprintf(fh,"\ttree gen.%d = [&U] %s\n", curstate->currentGeneration, curstate->tr->tree_string);

}



void printSample(state *curstate)
{
  exabayes_printTopology(curstate);
  printParams(curstate);
}



static void printSubsRates(state *prState ,int model, int numSubsRates)
{
  pInfo *partition = getPartition(prState, model); 

  assert(partition->dataType = DNA_DATA);
  int i;
  PRINT("Subs rates[%d]: ", model);
  for(i=0; i<numSubsRates; i++)
    PRINT("%d => %.3f, ", i, partition->substRates[i]);
    PRINT("\n");
  PRINT("frequencies[%d]: ", model);
  for(i=0; i<prState->frequRemem.numFrequRates; i++)
    PRINT("%d => %.3f, ", i, partition->frequencies[i]); 

  PRINT("\n");
}




void printIfPresent(proposalFunction *pf)
{
  int acc = pf->successCtr.acc,
    rejc = pf->successCtr.rej; 

  if(acc != 0 || rejc != 0)
    PRINT("%s: %d/%d (%.0f%%)\t", pf->name, acc , rejc,  ( (double)(acc) / (double)( (acc + rejc) + 0.0001))* 100   ); 
}



static void printHotChains(int runId)
{
  state *start =  gAInfo.allChains +  runId *  gAInfo.numberCoupledChains; 
  
  for(int i = 1; i < gAInfo.numberCoupledChains; ++i)
    {
      int index = 0; 
      for(int  j = 0; j < gAInfo.numberCoupledChains; ++j)
	{
	  state *chain =   start + j  ; 
	  if(chain->couplingId == i)
	    index = j; 
	}
      
      state *chain = start +index ;      

      double myHeat = getChainHeat(chain);
      assert(chain->couplingId < gAInfo.numberCoupledChains); 
      assert(chain->couplingId > 0 ) ; 
      assert( myHeat < 1.f);
      PRINT("lnl_beta(%.2f)=%.2f\t", myHeat, chain->likelihood); 

    }

  PRINT("successful swaps: %d", gAInfo.successFullSwitchesBatch); 
  gAInfo.successFullSwitchesBatch = 0; 

}


void chainInfo(state *chain)
{
  assert(chain->couplingId == 0) ; /* we are the cold chain   */

  int runId = chain->id / gAInfo.numberCoupledChains; 

  PRINT( "[run: %d] [TIME %.2f] gen: %d Likelihood: %.2f\t",runId,   gettime()  - timeIncrement  , chain->currentGeneration, chain->tr->likelihood);
  printHotChains(runId); 
  PRINT("\n"); 

  /* just output how much time has passed since the last increment */
  timeIncrement = gettime(); 	

  char* names[]= {"TOPO", "MODEL", "BL", "FREQ", "RATE_HET"} ;   

  for(int i = 1; i < NUM_PROP_CATS; ++i)
    {
      category_t type = i; 
      PRINT("%s:", names[i-1]);       
      for(int i = 0; i < chain->numProposals;++i )
	{
	  proposalFunction *pf = chain->proposals[i]; 	
	  if(type == pf->category)
	    {
	      PRINT("\t"); 
	      printIfPresent(pf);
	    }
	}
      PRINT("\n"); 
    }

  PRINT("\n");
}








#if 0 
/* QUARANTINE  */
void chainInfoOutput(state *curstate )
{

  PRINT( "[chain: %d] [TIME %.2f] gen: %d Likelihood: %f StartLH: %f \n",curstate->id,   gettime()  - timeIncrement  , curstate->currentGeneration, curstate->tr->likelihood, curstate->tr->startLH);

  /* just output how much time has passed since the last increment */
  timeIncrement = gettime(); 	
  
  PRINT( "Topo: eSpr: %d/%d (%d%%)\teSpr_mapped: %d/%d\t(%d%%) \n"
	 , curstate->acceptedProposals[E_SPR]	, curstate->rejectedProposals[E_SPR] , (int)(curstate->acceptedProposals[E_SPR]*100/
											     (curstate->acceptedProposals[E_SPR]+curstate->rejectedProposals[E_SPR]+0.0001)
),
	 curstate->acceptedProposals[E_SPR_MAPPED]	, curstate->rejectedProposals[E_SPR_MAPPED] , (int)(curstate->acceptedProposals[E_SPR_MAPPED]*100/(curstate->acceptedProposals[E_SPR_MAPPED]+curstate->rejectedProposals[E_SPR_MAPPED]+0.0001)));
  
  PRINT( "Model: slidingWindow: %d/%d (%d%%) biunif bin model: %d/%d (%d%%) biunif perm model: %d/%d (%d%%) single biunif model: %d/%d (%d%%) all biunif model: %d/%d (%d%%) \n",curstate->acceptedProposals[UPDATE_MODEL]	, curstate->rejectedProposals[UPDATE_MODEL] , (int)(curstate->acceptedProposals[UPDATE_MODEL]*100/(curstate->acceptedProposals[UPDATE_MODEL]+curstate->rejectedProposals[UPDATE_MODEL]+0.0001)) ,
	 curstate->acceptedProposals[UPDATE_MODEL_BIUNIF]	, curstate->rejectedProposals[UPDATE_MODEL_BIUNIF] , (int)(curstate->acceptedProposals[UPDATE_MODEL_BIUNIF]*100/(curstate->acceptedProposals[UPDATE_MODEL_BIUNIF]+curstate->rejectedProposals[UPDATE_MODEL_BIUNIF]+0.0001)) ,
	 curstate->acceptedProposals[UPDATE_MODEL_PERM_BIUNIF]	, curstate->rejectedProposals[UPDATE_MODEL_PERM_BIUNIF] , (int)(curstate->acceptedProposals[UPDATE_MODEL_PERM_BIUNIF]*100/(curstate->acceptedProposals[UPDATE_MODEL_PERM_BIUNIF]+curstate->rejectedProposals[UPDATE_MODEL_PERM_BIUNIF]+0.0001)) ,
	 curstate->acceptedProposals[UPDATE_MODEL_SINGLE_BIUNIF]	, curstate->rejectedProposals[UPDATE_MODEL_SINGLE_BIUNIF] , (int)(curstate->acceptedProposals[UPDATE_MODEL_SINGLE_BIUNIF]*100/(curstate->acceptedProposals[UPDATE_MODEL_SINGLE_BIUNIF]+curstate->rejectedProposals[UPDATE_MODEL_SINGLE_BIUNIF]+0.0001)) ,
	 curstate->acceptedProposals[UPDATE_MODEL_ALL_BIUNIF]	, curstate->rejectedProposals[UPDATE_MODEL_ALL_BIUNIF] , (int)(curstate->acceptedProposals[UPDATE_MODEL_ALL_BIUNIF]*100/(curstate->acceptedProposals[UPDATE_MODEL_ALL_BIUNIF]+curstate->rejectedProposals[UPDATE_MODEL_ALL_BIUNIF]+0.0001)));
  
  PRINT( "Frequencies: unif: %d/%d (%d%%) \n"
  , curstate->acceptedProposals[UPDATE_FREQUENCIES_BIUNIF]	, curstate->rejectedProposals[UPDATE_FREQUENCIES_BIUNIF] , (int)(curstate->acceptedProposals[UPDATE_FREQUENCIES_BIUNIF]*100/(curstate->acceptedProposals[UPDATE_FREQUENCIES_BIUNIF]+curstate->rejectedProposals[UPDATE_FREQUENCIES_BIUNIF]+0.0001)));
  
  PRINT( "Gamma: slidingWindow: %d/%d (%d%%) gaExp: %d/%d (%d%%) \n",curstate->acceptedProposals[UPDATE_GAMMA]	, curstate->rejectedProposals[UPDATE_GAMMA] , (int)(curstate->acceptedProposals[UPDATE_GAMMA]*100/(curstate->acceptedProposals[UPDATE_GAMMA]+curstate->rejectedProposals[UPDATE_GAMMA]+0.0001)),
	 curstate->acceptedProposals[UPDATE_GAMMA_EXP]	, curstate->rejectedProposals[UPDATE_GAMMA_EXP] , (int)(curstate->acceptedProposals[UPDATE_GAMMA_EXP]*100/(curstate->acceptedProposals[UPDATE_GAMMA_EXP]+curstate->rejectedProposals[UPDATE_GAMMA_EXP]+0.0001)));
  
  PRINT( "Branchlength: Slidingwindow: %d/%d (%d%%) blBiunif: %d/%d (%d%%) blExp: %d/%d (%d%%)\n",curstate->acceptedProposals[UPDATE_SINGLE_BL], curstate->rejectedProposals[UPDATE_SINGLE_BL], (int)(curstate->acceptedProposals[UPDATE_SINGLE_BL]*100/(curstate->acceptedProposals[UPDATE_SINGLE_BL]+curstate->rejectedProposals[UPDATE_SINGLE_BL]+0.0001)),
	 curstate->acceptedProposals[UPDATE_SINGLE_BL_EXP], curstate->rejectedProposals[UPDATE_SINGLE_BL_EXP], (int)(curstate->acceptedProposals[UPDATE_SINGLE_BL_EXP]*100/(curstate->acceptedProposals[UPDATE_SINGLE_BL_EXP]+curstate->rejectedProposals[UPDATE_SINGLE_BL_EXP]+0.0001)),
	 curstate->acceptedProposals[UPDATE_SINGLE_BL_BIUNIF], curstate->rejectedProposals[UPDATE_SINGLE_BL_BIUNIF], (int)(curstate->acceptedProposals[UPDATE_SINGLE_BL_BIUNIF]*100/(curstate->acceptedProposals[UPDATE_SINGLE_BL_BIUNIF]+curstate->rejectedProposals[UPDATE_SINGLE_BL_BIUNIF]+0.0001)));
  
  PRINT( "Hastings: %f new Prior: %f Current Prior: %f\n",curstate->hastings, curstate->newprior, curstate->curprior);
  
  if(getNumberOfPartitions(curstate->tr) < 10) 
    {
      for(int printModel=0; printModel< getNumberOfPartitions(curstate->tr);printModel++) 
	printSubsRates(curstate, printModel, curstate->modelRemem.numSubsRates);
    }
    
  PRINT("\n");
  //printSubsRates(curstate, curstate->modelRemem.model, curstate->modelRemem.numSubsRates);
//printf("numSubs: %d numFrequ: %d\n\n",curstate->modelRemem.numSubsRates,curstate->frequRemem.numFrequRates);

}
#endif

