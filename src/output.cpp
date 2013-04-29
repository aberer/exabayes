#include <sstream>

#include "axml.h"
#include "bayes.h"
#include "main-common.h"
#include "globals.h"
#include "output.h"
#include "adapters.h"
#include "topology-utils.h"
#include "TreeAln.hpp"


// #include "chain.h"

#include "Chain.hpp"




/**
   @brief prints the tree as a debug message. 
*/
void debug_printTree(Chain *chain)
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
void debug_printNodeEnvironment(Chain *chain, int nodeID )
{
#ifdef DEBUG_SHOW_TREE
  nodeptr 
    p =  chain->tr->nodep[nodeID]; 

  printInfo(chain, "STATE node %d:\thooked to %d (bl=%f),\t%d (bl=%f)\tand %d (bl=%f)\n", p->number, p->back->number, traln->getBranchLength( p->back->number,0), 
	    p->next->back->number, traln->getBranchLength( p->next->back->number,0),
	    p->next->next->back->number, traln->getBranchLength( p->next->next->back->number,0)
	    );   
#endif
}



void debug_printAccRejc(Chain *chain, proposalFunction *pf, boolean accepted) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  if(isOutputProcess())
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
      if(isOutputProcess())
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
void printInfo(Chain *chain, const char *format, ...)
{  
  if( isOutputProcess())
    {
      printf("[run %d / heat %d / gen %d] ", chain->id / gAInfo.numberCoupledChains, chain->couplingId, chain->currentGeneration); 
      va_list args;
      va_start(args, format);     
      vprintf(format, args );
      va_end(args);
    }
}


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

  sprintf(treestr, ":%.7f", branchLengthToReal(tr, p->z[0]));
  
  while (*treestr) treestr++;
  return  treestr;
}


/* TODO maybe generate random run id */


static void printNexusTreeFileStart(Chain *chain)
{
  FILE *fh = chain->topologyFile; 

  fprintf(fh, "#NEXUS\n[ID: %s]\n[Param: tree]\nbegin trees;\n\ttranslate\n", "TODO" ); 

  for(int i = 0; i < chain->traln->getTr()->mxtips-1; ++i)
    fprintf(fh, "\t\t%d %s,\n", i+1, chain->traln->getTr()->nameList[i+1]); 
  fprintf(fh, "\t\t%d %s;\n", chain->traln->getTr()->mxtips, chain->traln->getTr()->nameList[chain->traln->getTr()->mxtips]);
}


static void printParamFileStart(Chain *chain)
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

static void printParams(Chain *chain)
{
  tree *tr = chain->traln->getTr(); 

  FILE *fh = chain->outputParamFile; 

  double treeLength = branchLengthToReal(tr, getTreeLength(chain->traln,tr->nodep[1]->back )); 
  assert(treeLength != 0.); 
  fprintf(fh, "%d\t%f\t%.3f", chain->currentGeneration,
	  tr->likelihood,  
	  treeLength); 

  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      pInfo *partition = chain->traln->getPartition(i); 
      fprintf(fh, "\t%f\t%f", chain->traln->accessPartitionLH( i),chain->traln->getAlpha(i)) ; 
      for(int j = 0; j < 6 ; ++j) /* TODO */
	fprintf(fh, "\t%.2f", partition->substRates[j]); 
      for(int j = 0; j < 4 ; ++j) /* TODO */
	fprintf(fh, "\t%.2f", chain->traln->getFrequency(i,j));
    }

  fprintf(fh, "\n"); 
  fflush(fh); 
}



void initializeOutputFiles(Chain *chain)
{
  printNexusTreeFileStart(chain); 
  printParamFileStart(chain);
}


void finalizeOutputFiles(Chain *chain)
{
  FILE *fh = chain->topologyFile; 

  /* topo file  */
  fprintf(fh, "end;\n"); 
  fclose(fh); 
}



  /* TODO what about per model brach lengths? how does mrB do this? */
static void exabayes_printTopology(Chain *chain)
{  
  assert(chain->couplingId == 0);
  tree *tr = chain->traln->getTr();
  FILE *fh = chain->topologyFile; 
  memset(chain->traln->getTr()->tree_string, 0, chain->traln->getTr()->treeStringLength * sizeof(char) ); 
  
  Tree2stringNexus(chain->traln->getTr()->tree_string, tr,  chain->traln->getTr()->start->back, 0 ); 
  fprintf(fh,"\ttree gen.%d = [&U] %s\n", chain->currentGeneration, chain->traln->getTr()->tree_string);

}



void printSample(Chain *chain)
{
  exabayes_printTopology(chain);
  printParams(chain);
}




void printIfPresent(proposalFunction *pf)
{
  stringstream buf; 
  buf << pf->sCtr; 
  PRINT("%s: %s\t", pf->name, buf.str().c_str());   
}

static void printHotChains(int runId)
{
  Chain *start =  gAInfo.allChains +  runId *  gAInfo.numberCoupledChains; 

  for(int i = 1; i < gAInfo.numberCoupledChains; ++i)
    {
      int index = 0; 
      for(int  j = 0; j < gAInfo.numberCoupledChains; ++j)
	{
	  Chain *chain =   start + j  ; 
	  if(chain->couplingId == i)
	    index = j; 
	}
      
      Chain *chain = start +index ;      

      double myHeat = chain->getChainHeat();

      assert(chain->couplingId < gAInfo.numberCoupledChains); 
      assert(chain->couplingId > 0 ) ; 
      assert( myHeat < 1.f);

      PRINT("lnl_beta(%.2f)=%.2f\t", myHeat, chain->traln->getTr()->likelihood); 
    }
}



/**
   @brief dumps infos on the Chain of the chain

   @param chain -- the pointer to the first chain that belongs to a particular run in the array of chains 
 */
void chainInfo(Chain *chain)
{
  assert(chain->couplingId == 0) ; /* we are the cold chain   */

  int runId = chain->id / gAInfo.numberCoupledChains; 
  tree *tr = chain->traln->getTr(); 

  PRINT( "[run: %d] [time %.2f] gen: %d Likelihood: %.2f\tTL=%.2f\t",runId,   gettime()  - timeIncrement  , chain->currentGeneration, chain->traln->getTr()->likelihood, branchLengthToReal(tr, getTreeLength(chain->traln, tr->nodep[1]->back)));
  printHotChains(runId); 
  // printSwapInfo(runId);   
  // TODO 
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

