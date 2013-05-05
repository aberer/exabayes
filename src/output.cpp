#include <sstream>

#include "axml.h"
#include "bayes.h"
// #include "main-common.h"
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
void debug_printTree(TreeAln &traln)
{
#ifdef DEBUG_SHOW_TREE  
  tree *tr = traln.getTr();
  Tree2stringNexus(tr->tree_string, tr, tr->start->back, 0); 
  cout << "tree: "<< tr->tree_string << endl; 
#endif
}




/**
   @brief prints the entire environment of a node with branch lengths
 */
void debug_printNodeEnvironment(TreeAln &traln, int nodeID )
{
#ifdef DEBUG_SHOW_TREE
  tree *tr = traln.getTr();
  nodeptr 
    p =  tr->nodep[nodeID]; 

  
  printf("STATE node %d:\thooked to %d (bl=%f),\t%d (bl=%f)\tand %d (bl=%f)\n", p->number, p->back->number, traln.getBranchLength( p->back,0), 
	 p->next->back->number, traln.getBranchLength( p->next->back,0),
	    p->next->next->back->number, traln.getBranchLength( p->next->next->back,0)
	    );   
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
  fflush(fh);
}



void printSample(Chain *chain)
{
  exabayes_printTopology(chain);
  printParams(chain);
}



