#include "globals.h"
#include "common.h"


#include "config.h"
#include "axml.h"
#include "proposals.h"
#include "main-common.h"


#include "output.h"



extern double masterTime; 

void makeFileNames(void)
{
  int infoFileExists = 0;

  strcpy(topologyFile, workdir); 
  strcat(topologyFile, PROGRAM_NAME);
  strcat(topologyFile, "_topologies."); 
  strcat(topologyFile, run_id); 
 
  strcpy(outputParamFile,workdir); 
  strcat(outputParamFile, PROGRAM_NAME);
  strcat(outputParamFile, "_parameters."); 
  strcat(outputParamFile, run_id); 
  
  strcpy(infoFileName, workdir); 
  strcat(infoFileName, PROGRAM_NAME);
  strcat(infoFileName, "_infoFile."); 
  strcat(infoFileName, run_id); 


  strcpy(binaryChainState, workdir); 
  strcat(binaryChainState, PROGRAM_NAME);
  strcat(binaryChainState, "_binChainState."); 
  strcat(binaryChainState, run_id); 

  infoFileExists = filexists(infoFileName);


  if(infoFileExists)
  {
    printf("RAxML output files with the run ID <%s> already exist \n", run_id);
    printf("in directory %s ...... exiting\n", workdir);
#ifdef PRODUCTIVE
    exit(-1);
#endif
  }
}



double exabayes_getBranchLength(tree *tr, int perGene, nodeptr p)
{
  double 
    z = 0.0,
    x = 0.0;

  assert(perGene != NO_BRANCHES);
	      
  if(tr->numBranches == 1)
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
		      
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      assert(tr->partitionContributions[i] != -1.0);
	      assert(tr->fracchanges[i] != -1.0);
	      z = p->z[i];
	      if(z < zmin) 
		z = zmin;      	 
	      x = -log(z) * tr->fracchanges[i];
	      avgX += x * tr->partitionContributions[i];
	    }

	  x = avgX;
	}
      else
	{	
	  assert(tr->fracchanges[perGene] != -1.0);
	  assert(perGene >= 0 && perGene < tr->numBranches);
	  
	  z = p->z[perGene];
	  
	  if(z < zmin) 
	    z = zmin;      	 
	  
	  x = -log(z) * tr->fracchanges[perGene];	  
	}
    }

  return x;
}


static char *Tree2stringNexus(char *treestr, tree *tr, nodeptr p, int perGene )
{      
  if(isTip(p->number, tr->mxtips)) 
    {	       	        
      sprintf(treestr, "%d", p->number); 
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = Tree2stringNexus(treestr, tr, p->next->back, perGene);
      *treestr++ = ',';
      treestr = Tree2stringNexus(treestr, tr, p->next->next->back, perGene);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2stringNexus(treestr, tr, p->back, perGene);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back) 
    {
      sprintf(treestr, ";"); 
      return treestr;
    }

  sprintf(treestr, ":%8.20f", exabayes_getBranchLength(tr, perGene, p));
  
  while (*treestr) treestr++;
  return  treestr;
}


/* TODO maybe generate random run id */


static void printNexusTreeFileStart(state *curstate)
{
  FILE *fh = myfopen(topologyFile, "w"); 
  fprintf(fh, "#NEXUS\n[ID: %s]\n[Param: tree]\nbegin trees;\n\ttranslate\n", "TODO" ); 

  for(int i = 0; i < curstate->tr->mxtips-1; ++i)
    fprintf(fh, "\t\t%d %s,\n", i+1, curstate->tr->nameList[i+1]); 
  fprintf(fh, "\t\t%d %s;\n", curstate->tr->mxtips, curstate->tr->nameList[curstate->tr->mxtips]);
  fclose(fh); 
}


static void printParamFileStart(state *curstate)
{
  FILE *fh = myfopen(outputParamFile, "w"); 
  
  char *tmp = "TODO"; 
  fprintf(fh, "[ID: %s]\n", tmp);
  fprintf(fh, "Gen\tLnL\tTL\
\tr(A<->C)\tr(A<->G)\tr(A<->T)\tr(C<->G)\tr(C<->T)\tr(G<->T)\t\
pi(A)\tpi(C)\tpi(G)\tpi(T)\
\talpha\tpinvar\n"); 

  fclose(fh); 
}

static void printParams(state *curstate)
{
  int model = 0; 
  
  FILE *fh = myfopen(outputParamFile, "a"); 

  /* TODO tree length */
  fprintf(fh, "%d\t%f\t%f", curstate->currentGeneration,curstate->tr->likelihood, -1. ); 

  /* TODO what about multiple models?  */
  /* TODO that will not work in the future ...  */
  for(int i = 0; i< curstate->modelRemem.numSubsRates;  ++i )
    fprintf(fh, "\t%f" , curstate->tr->partitionData[model].substRates[i]); 

  
  /* TODO works on dna only currently  */
  for(int i = 0; i < 4; ++i)
    fprintf(fh, "\t%f", curstate->tr->partitionData[model].frequencies[i]); 

  /* TODO what is pinvar?  */
  fprintf(fh, "\t%f\t%f", curstate->tr->partitionData[model].alpha, -1.  ); 

  fprintf(fh, "\n");

  fclose(fh); 
}



void initializeOutputFiles(state *curstate)
{
  printNexusTreeFileStart(curstate); 
  printParamFileStart(curstate);
}


void finalizeOutputFiles(state *curstate)
{
  /* topo file  */
  FILE *fh = myfopen(topologyFile, "a"); 
  fprintf(fh, "end;\n"); 
  fclose(fh); 
  
  /* info file  */
PRINT("\n\nTotal execution time for %d generations: %f\n",  curstate->numGen, gettime() - masterTime); 
}



  /* TODO what about per model brach lengths? how does mrB do this? */
static void exabayes_printTopology(state *curstate)
{
  FILE *fh = myfopen(topologyFile, "a");  
  memset(curstate->tr->tree_string, 0, curstate->tr->treeStringLength * sizeof(char) ); 
  
  Tree2stringNexus(curstate->tr->tree_string, curstate->tr,  curstate->tr->start->back, 0 ); 
  fprintf(fh,"\ttree gen.%d = [&U] %s\n", curstate->currentGeneration, curstate->tr->tree_string);

  fclose(fh);   
}



void printSample(state *curstate)
{
  exabayes_printTopology(curstate);
  printParams(curstate);
}



static void printSubsRates(tree *tr,int model, int numSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  PRINT("Subs rates: ");
  for(i=0; i<numSubsRates; i++)
    PRINT("%d => %.3f, ", i, tr->partitionData[model].substRates[i]);
  PRINT("\n\n");
}


/* TODO modify this */
void chainInfoOutput(state *curstate, int sum_radius_accept, int sum_radius_reject )
{
  PRINT( "propb: %d %f %f spr: %d (%d) model: %d (%d) ga: %d (%d) bl: %d (%d) blExp: %d (%d) %f %f %f radius: %f %f\n",
	 curstate->currentGeneration, curstate->tr->likelihood, curstate->tr->startLH, 
	 curstate->acceptedProposals[SPR]	, curstate->rejectedProposals[SPR] ,
	 curstate->acceptedProposals[UPDATE_MODEL]	, curstate->rejectedProposals[UPDATE_MODEL] ,
	 curstate->acceptedProposals[UPDATE_GAMMA]	, curstate->rejectedProposals[UPDATE_GAMMA] ,
	 curstate->acceptedProposals[UPDATE_SINGLE_BL], curstate->rejectedProposals[UPDATE_SINGLE_BL],
	 curstate->acceptedProposals[UPDATE_SINGLE_BL_EXP], curstate->rejectedProposals[UPDATE_SINGLE_BL_EXP],
	 curstate->hastings, curstate->newprior, curstate->curprior,
	 sum_radius_accept / (float)curstate->acceptedProposals[SPR],
	 sum_radius_reject / (float)curstate->rejectedProposals[SPR] );
  
  printSubsRates(curstate->tr, curstate->modelRemem.model, curstate->modelRemem.numSubsRates);

#ifdef WITH_PERFORMANCE_MEASUREMENTS
  perf_timer_print( &all_timer );
#endif     
}



/* TODO */
void newChainInfoOutput(state *curstate)
{
  /* current values */

  /* acceptance  / rejection values  */ 
}

