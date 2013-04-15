


#include "axml.h"
#include "bayes.h"
#include "main-common.h"
#include "globals.h"
#include "chain.h"
#include "output.h"
#include "adapters.h"
#include "topology-utils.h"
#include "TreeAln.hpp"



/**
   @brief prints the tree as a debug message. 
*/
void debug_printTree(state *chain)
{
#ifdef DEBUG_SHOW_TREE  
  char tmp[100000];
  Tree2stringNexus(tmp, chain, chain->tr->start->back, 0); 
  printInfo(chain, "tree: %s\n", tmp); 
#endif
}




/**
   @brief prints the entire environment of a node with branch lengths
 */
void debug_printNodeEnvironment(state *chain, int nodeID )
{
#ifdef DEBUG_SHOW_TREE
  nodeptr 
    p =  chain->tr->nodep[nodeID]; 

  printInfo(chain, "STATE node %d:\thooked to %d (bl=%f),\t%d (bl=%f)\tand %d (bl=%f)\n", p->number, p->back->number, p->back->z[0], 
	    p->next->back->number, p->next->back->z[0],
	    p->next->next->back->number, p->next->next->back->z[0]
	    );   
#endif
}



void debug_printAccRejc(state *chain, proposalFunction *pf, boolean accepted) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  if(processID == 0)
    {
      if(accepted)
	printInfo(chain, "ACCEPTING\t");   
      else 
	printInfo(chain, "rejecting\t");   	  
      printf("%s %g\n" ,pf->name, chain->traln->getTr()->likelihood); 
    }
#endif
}


void debug_checkTreeConsistency(tree *tr)
{
#ifdef DEBUG_CHECK_TREE_CONSISTENCY
  // tree *tr = chain->tr; 
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
}

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




void printOrientation(tree *tr, nodeptr p)
{
  if(isTip(p->number,tr->mxtips))
    printf("%d is tip\n", p->number); 
  else if(p->x )
    printf("%d is orientated towards %d\n", p->number,p->back->number); 
  else if(p->next->x)
    printf("%d is orientated towards %d\n", p->number,p->next->back->number); 
  else if(p->next->next->x)
    printf("%d is orientated towards %d\n", p->number,p->next->next->back->number); 
  else 
    assert(0); 
}




/**
   @brief prints a message with associated chain/run/heat information. 

   Please always use this, when possible, it also assures that the
   message is only printed once.
 */ 
void printInfo(state *chain, const char *format, ...)
{  
  if(
#if HAVE_PLL == 1
     true || 
#endif
     processID == 0)
    {
      printf("[run %d / heat %d / gen %d] ", chain->id / gAInfo.numberCoupledChains, chain->couplingId, chain->currentGeneration); 
      va_list args;
      va_start(args, format);     
      vprintf(format, args );
      va_end(args);
    }
}


/* double exabayes_getBranchLength(state *chain, int perGene, nodeptr p) */
/* { */
/*   tree *tr = chain->tr;  */

/*   double  */
/*     z = 0.0, */
/*     x = 0.0; */

/*   assert(perGene != NO_BRANCHES); */
	      
/*   if(NOT hasPergeneBL(tr)) */
/*     { */
/*       assert(tr->fracchange != -1.0); */
/*       z = p->z[0]; */
/*       if (z < zmin)  */
/* 	z = zmin;      	  */
      
/*       x = -log(z) * tr->fracchange;  */
/*     } */
/*   else */
/*     { */
/*       if(perGene == SUMMARIZE_LH) */
/* 	{ */
/* 	  int  */
/* 	    i; */
	  
/* 	  double  */
/* 	    avgX = 0.0; */
		      
/* 	  for(i = 0; i < getNumBranches(tr); i++) */
/* 	    { */

/* 	      assert(getPcontr(chain,i) != -1.0); */
/* 	      assert(getFracChange(chain,i) != -1.0); */

/* 	      z = p->z[i]; */
/* 	      if(z < zmin)  */
/* 		z = zmin;      	  */
/* 	      x = -log(z) * getFracChange(chain,i); */
/* 	      avgX += x * getPcontr(chain,i); */
/* 	    } */

/* 	  x = avgX; */
/* 	} */
/*       else */
/* 	{	 */
/* 	  assert(getFracChange(chain,perGene) != -1.0); */
/* 	  assert(perGene >= 0 && perGene < getNumBranches(tr)); */
	  
/* 	  z = p->z[perGene]; */
	  
/* 	  if(z < zmin)  */
/* 	    z = zmin;      	  */
	  
/* 	  x = -log(z) * getFracChange(chain,perGene); */
/* 	} */
/*     } */

/*   return x; */
/* } */


char *Tree2stringNexus(char *treestr, tree *tr , nodeptr p, int perGene )
{      
  // tree *tr = chain->tr; 

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

  sprintf(treestr, ":%g", branchLengthToReal(tr, p->z[0]));
  
  while (*treestr) treestr++;
  return  treestr;
}


/* TODO maybe generate random run id */


static void printNexusTreeFileStart(state *chain)
{
  FILE *fh = chain->topologyFile; 

  fprintf(fh, "#NEXUS\n[ID: %s]\n[Param: tree]\nbegin trees;\n\ttranslate\n", "TODO" ); 

  for(int i = 0; i < chain->traln->getTr()->mxtips-1; ++i)
    fprintf(fh, "\t\t%d %s,\n", i+1, chain->traln->getTr()->nameList[i+1]); 
  fprintf(fh, "\t\t%d %s;\n", chain->traln->getTr()->mxtips, chain->traln->getTr()->nameList[chain->traln->getTr()->mxtips]);
}


static void printParamFileStart(state *chain)
{
  FILE *fh = chain->outputParamFile; 

  char *tmp = "TODO"; 
  fprintf(fh, "[ID: %s]\n", tmp);
  fprintf(fh, "Gen\tLnL\tTL"); 

  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      fprintf(fh, "\tlnl.%d\talhpa.%d", i,i); 
      fprintf(fh, "\tr.%d(A<->C)\tr.%d(A<->G)\tr.%d(A<->T)\tr.%d(C<->G)\tr.%d(C<->T)\tr.%d(G<->T)",i, i,i,i,i,i); 
      fprintf(fh, "\tpi(A).%d\tpi(C).%d\tpi(G).%d\tpi(T).%d", i,i,i,i); 
    }

  fprintf(fh, "\n"); 
  fflush(fh); 
}

static void printParams(state *chain)
{
  tree *tr = chain->traln->getTr(); 

  FILE *fh = chain->outputParamFile; 

  double treeLength = branchLengthToReal(tr, getTreeLength(tr,tr->nodep[1]->back )); 
  assert(treeLength != 0.); 
  fprintf(fh, "%d\t%f\t%.3f", chain->currentGeneration,
	  tr->likelihood,  
	  treeLength); 

  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      pInfo *partition = chain->traln->getPartition(i); 
      fprintf(fh, "\t%f\t%f", chain->traln->accessPartitionLH( i),partition->alpha) ; 
      for(int j = 0; j < 6 ; ++j) /* TODO */
	fprintf(fh, "\t%.2f", partition->substRates[j]); 
      for(int j = 0; j < 4 ; ++j) /* TODO */
	fprintf(fh, "\t%.2f", partition->frequencies[j]) ;
    }

  fprintf(fh, "\n"); 
  fflush(fh); 
}



void initializeOutputFiles(state *chain)
{
  printNexusTreeFileStart(chain); 
  printParamFileStart(chain);
}


void finalizeOutputFiles(state *chain)
{
  FILE *fh = chain->topologyFile; 

  /* topo file  */
  fprintf(fh, "end;\n"); 
  fclose(fh); 
}



  /* TODO what about per model brach lengths? how does mrB do this? */
static void exabayes_printTopology(state *chain)
{
  tree *tr = chain->traln->getTr();
  /* FILE *fh = myfopen(topologyFile, "a");   */
  FILE *fh = chain->topologyFile; 
  memset(chain->traln->getTr()->tree_string, 0, chain->traln->getTr()->treeStringLength * sizeof(char) ); 
  
  Tree2stringNexus(chain->traln->getTr()->tree_string, tr,  chain->traln->getTr()->start->back, 0 ); 
  fprintf(fh,"\ttree gen.%d = [&U] %s\n", chain->currentGeneration, chain->traln->getTr()->tree_string);

}



void printSample(state *chain)
{
  exabayes_printTopology(chain);
  printParams(chain);
}




#if 0 
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
#endif




void printIfPresent(proposalFunction *pf)
{
  double ratio = getRatioOverall(&(pf->sCtr)); 
  int acc = pf->sCtr.gAcc; 
  int rej = pf->sCtr.gRej; 

  if(acc != 0 || rej != 0)
    PRINT("%s: %d/%d (%.1f%%)\t", pf->name, acc  , rej ,  ratio * 100     ); 
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

      PRINT("lnl_beta(%.2f)=%.2f\t", myHeat, chain->traln->getTr()->likelihood); 
    }
}




static void printSwapInfo(int runId)
{
  successCtr
    *ctrMatrix = gAInfo.swapInfo[runId]; 
  
  int cnt = 0; 
  for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
    {
      if(i < gAInfo.numberCoupledChains - 1 )
	PRINT("("); 

      for(int j = 0; j < gAInfo.numberCoupledChains; ++j)
	{
	  successCtr *ctr = & ( ctrMatrix[cnt]) ; 
	  if(i < j )
	    {	      
	      PRINT("%.1f%%,", i,j, 100 * getRatioOverall(ctr));
	    }
	  else 
	    {
	      assert( ctr->gAcc == 0 
			 && ctr->gRej == 0
			 && ctr->lAcc == 0
			 && ctr->lRej == 0);  
	    }
	  cnt++; 
	}
      if(i < gAInfo.numberCoupledChains - 1 )
	PRINT(")"); 
    }

  /* PRINT("\tcurrentRatio(0,1)=%.1f%%", getRatioLocal( ctrMatrix + 1) );  */
}




/**
   @brief dumps infos on the state of the chain

   @param chain -- the pointer to the first chain that belongs to a particular run in the array of chains 
 */
void chainInfo(state *chain)
{
  assert(chain->couplingId == 0) ; /* we are the cold chain   */

  int runId = chain->id / gAInfo.numberCoupledChains; 
  tree *tr = chain->traln->getTr(); 

  PRINT( "[run: %d] [time %.2f] gen: %d Likelihood: %.2f\tTL=%.2f\t",runId,   gettime()  - timeIncrement  , chain->currentGeneration, chain->traln->getTr()->likelihood, branchLengthToReal(tr, getTreeLength(tr, tr->nodep[1]->back)));
  printHotChains(runId); 
  printSwapInfo(runId);   
  PRINT("\n"); 

  /* just output how much time has passed since the last increment */
  timeIncrement = gettime(); 	

  char* names[]= {"TOPO", "BL", "FREQ", "MODEL", "HETER"} ;   

  for(int i = 1; i < NUM_PROP_CATS + 1 ; ++i)
    {
      category_t type = category_t(i); 
      PRINT("%s:", names[i-1]);       
      for(int j = 0; j < chain->numProposals;++j )
	{
	  proposalFunction *pf = chain->proposals[j]; 	
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
void chainInfoOutput(state *chain )
{

  PRINT( "[chain: %d] [TIME %.2f] gen: %d Likelihood: %f StartLH: %f \n",chain->id,   gettime()  - timeIncrement  , chain->currentGeneration, chain->traln->getTr()->likelihood, chain->traln->getTr()->startLH);

  /* just output how much time has passed since the last increment */
  timeIncrement = gettime(); 	
  
  PRINT( "Topo: eSpr: %d/%d (%d%%)\teSpr_mapped: %d/%d\t(%d%%) \n"
	 , chain->acceptedProposals[E_SPR]	, chain->rejectedProposals[E_SPR] , (int)(chain->acceptedProposals[E_SPR]*100/
											     (chain->acceptedProposals[E_SPR]+chain->rejectedProposals[E_SPR]+0.0001)
),
	 chain->acceptedProposals[E_SPR]	, chain->rejectedProposals[E_SPR_MAPPED] , (int)(chain->acceptedProposals[E_SPR_MAPPED]*100/(chain->acceptedProposals[E_SPR_MAPPED]+chain->rejectedProposals[E_SPR_MAPPED]+0.0001)));
  
  PRINT( "Model: slidingWindow: %d/%d (%d%%) biunif bin model: %d/%d (%d%%) biunif perm model: %d/%d (%d%%) single biunif model: %d/%d (%d%%) all biunif model: %d/%d (%d%%) \n",chain->acceptedProposals[UPDATE_MODEL]	, chain->rejectedProposals[UPDATE_MODEL] , (int)(chain->acceptedProposals[UPDATE_MODEL]*100/(chain->acceptedProposals[UPDATE_MODEL]+chain->rejectedProposals[UPDATE_MODEL]+0.0001)) ,
	 chain->acceptedProposals[UPDATE_MODEL_BIUNIF]	, chain->rejectedProposals[UPDATE_MODEL_BIUNIF] , (int)(chain->acceptedProposals[UPDATE_MODEL_BIUNIF]*100/(chain->acceptedProposals[UPDATE_MODEL_BIUNIF]+chain->rejectedProposals[UPDATE_MODEL_BIUNIF]+0.0001)) ,
	 chain->acceptedProposals[UPDATE_MODEL_PERM_BIUNIF]	, chain->rejectedProposals[UPDATE_MODEL_PERM_BIUNIF] , (int)(chain->acceptedProposals[UPDATE_MODEL_PERM_BIUNIF]*100/(chain->acceptedProposals[UPDATE_MODEL_PERM_BIUNIF]+chain->rejectedProposals[UPDATE_MODEL_PERM_BIUNIF]+0.0001)) ,
	 chain->acceptedProposals[UPDATE_MODEL_SINGLE_BIUNIF]	, chain->rejectedProposals[UPDATE_MODEL_SINGLE_BIUNIF] , (int)(chain->acceptedProposals[UPDATE_MODEL_SINGLE_BIUNIF]*100/(chain->acceptedProposals[UPDATE_MODEL_SINGLE_BIUNIF]+chain->rejectedProposals[UPDATE_MODEL_SINGLE_BIUNIF]+0.0001)) ,
	 chain->acceptedProposals[UPDATE_MODEL_ALL_BIUNIF]	, chain->rejectedProposals[UPDATE_MODEL_ALL_BIUNIF] , (int)(chain->acceptedProposals[UPDATE_MODEL_ALL_BIUNIF]*100/(chain->acceptedProposals[UPDATE_MODEL_ALL_BIUNIF]+chain->rejectedProposals[UPDATE_MODEL_ALL_BIUNIF]+0.0001)));
  
  PRINT( "Frequencies: unif: %d/%d (%d%%) \n"
  , chain->acceptedProposals[UPDATE_FREQUENCIES_BIUNIF]	, chain->rejectedProposals[UPDATE_FREQUENCIES_BIUNIF] , (int)(chain->acceptedProposals[UPDATE_FREQUENCIES_BIUNIF]*100/(chain->acceptedProposals[UPDATE_FREQUENCIES_BIUNIF]+chain->rejectedProposals[UPDATE_FREQUENCIES_BIUNIF]+0.0001)));
  
  PRINT( "Gamma: slidingWindow: %d/%d (%d%%) gaExp: %d/%d (%d%%) \n",chain->acceptedProposals[UPDATE_GAMMA]	, chain->rejectedProposals[UPDATE_GAMMA] , (int)(chain->acceptedProposals[UPDATE_GAMMA]*100/(chain->acceptedProposals[UPDATE_GAMMA]+chain->rejectedProposals[UPDATE_GAMMA]+0.0001)),
	 chain->acceptedProposals[UPDATE_GAMMA_EXP]	, chain->rejectedProposals[UPDATE_GAMMA_EXP] , (int)(chain->acceptedProposals[UPDATE_GAMMA_EXP]*100/(chain->acceptedProposals[UPDATE_GAMMA_EXP]+chain->rejectedProposals[UPDATE_GAMMA_EXP]+0.0001)));
  
  PRINT( "Branchlength: Slidingwindow: %d/%d (%d%%) blBiunif: %d/%d (%d%%) blExp: %d/%d (%d%%)\n",chain->acceptedProposals[UPDATE_SINGLE_BL], chain->rejectedProposals[UPDATE_SINGLE_BL], (int)(chain->acceptedProposals[UPDATE_SINGLE_BL]*100/(chain->acceptedProposals[UPDATE_SINGLE_BL]+chain->rejectedProposals[UPDATE_SINGLE_BL]+0.0001)),
	 chain->acceptedProposals[UPDATE_SINGLE_BL_EXP], chain->rejectedProposals[UPDATE_SINGLE_BL_EXP], (int)(chain->acceptedProposals[UPDATE_SINGLE_BL_EXP]*100/(chain->acceptedProposals[UPDATE_SINGLE_BL_EXP]+chain->rejectedProposals[UPDATE_SINGLE_BL_EXP]+0.0001)),
	 chain->acceptedProposals[UPDATE_SINGLE_BL_BIUNIF], chain->rejectedProposals[UPDATE_SINGLE_BL_BIUNIF], (int)(chain->acceptedProposals[UPDATE_SINGLE_BL_BIUNIF]*100/(chain->acceptedProposals[UPDATE_SINGLE_BL_BIUNIF]+chain->rejectedProposals[UPDATE_SINGLE_BL_BIUNIF]+0.0001)));
  
  PRINT( "Hastings: %f new Prior: %f Current Prior: %f\n",chain->hastings, chain->newprior, chain->curprior);
  
  if(getNumberOfPartitions(chain->traln->getTr()) < 10) 
    {
      for(int printModel=0; printModel< getNumberOfPartitions(chain->traln->getTr());printModel++) 
	printSubsRates(chain, printModel, chain->modelRemem.numSubsRates);
    }
    
  PRINT("\n");
  //printSubsRates(chain, chain->modelRemem.model, chain->modelRemem.numSubsRates);
//printf("numSubs: %d numFrequ: %d\n\n",chain->modelRemem.numSubsRates,chain->frequRemem.numFrequRates);

}
#endif

