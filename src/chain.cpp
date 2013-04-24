
#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "main-common.h"
#include "chain.h"
#include "proposals.h"
#include "output.h"
#include "treeRead.h"
#include "topology-utils.h"
#include "randomness.h"
#include "eval.h"
#include "adapters.h"
#include "randomTree.h" 

#include "TreeAln.hpp"
#include "LnlRestorer.hpp"
#include "Topology.hpp"

#include "BipartitionHash.hpp"


#include "CoupledChains.hpp"


#include <vector>
using namespace std; 


// void makeChainFileNames(Chain *chain, int num)
// {

// }


/**
   @brief corrects branch lengths when those are read from a file
*/
void traverseInitCorrect(nodeptr p, int *count, TreeAln *traln )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;

  for( i = 0; i < traln->getNumBranches(); i++)
    {
      double val = traln->getBranchLength(p, i); 
      traln->setBranchLengthSave(exp( - val  / tr->fracchange), i,p); 
    }
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitCorrect(q->back, count, traln);
	  q = q->next;
	} 
    }
}

void initParamDump(TreeAln *traln, paramDump *dmp)
{  
  tree *tr = traln->getTr();
  
  dmp->topology = new Topology(tr->mxtips); 
  dmp->infoPerPart = (perPartitionInfo*)exa_calloc(traln->getNumberOfPartitions(), sizeof(perPartitionInfo));
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *p =  dmp->infoPerPart + i ; 
      p->alpha = 0.5; 
      for(int j = 0; j < 6; ++j)
	p->substRates[j] = 0.5; 
      p->substRates[5] = 1; 
      for(int j = 0; j < 4; ++j)
	p->frequencies[j] = 0.25;       
    }

  dmp->branchLengths = (double*)exa_calloc(2 * tr->mxtips, sizeof(double));
  for(int i = 0; i < 2 * tr->mxtips; ++i)
    dmp->branchLengths[i] = INIT_BRANCH_LENGTHS;   
}


void copyParamDump(TreeAln *traln,paramDump *dest, const paramDump *src)
{
  /* no topo or branch lengths! this should be done by copytopo */
  
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *pDest =  dest->infoPerPart + i ,
	*pSrc = src->infoPerPart + i; 
      
      pDest->alpha = pSrc->alpha; 
      
      for(int j = 0; j < 6; ++j)
	pDest->substRates[j] = pSrc->substRates[j]; 
      for(int j = 0; j < 4; ++j)
	pDest->frequencies[j] = pSrc->frequencies[j]; 
    }

}




void copyState(Chain *dest, const Chain *src )
{
  TreeAln &thisTraln =  *(dest->traln); 
  TreeAln &rhsTraln = *(src->traln); 

  thisTraln = rhsTraln; 
  *(dest->traln ) = *(src->traln); 
  
  dest->dump.topology->saveTopology(thisTraln);
}


/**
   @brief An overloaded initialization function for all the chains. 

   This function should initialize all chain data structures and parse
   the config file.

   This function also decides which aln,tr structures are assigned to
   which chains.
 */ 
void initializeIndependentChains( analdef *adef, int seed, vector<CoupledChains*> &runs, initParamStruct *initParams )
{
  FILE *treeFH = NULL; 
  if( gAInfo.numberOfStartingTrees > 0 )
    treeFH = myfopen(tree_file, "r"); 

  int totalNumChains = gAInfo.numberOfRuns * gAInfo.numberCoupledChains;   
  PRINT("number of independent runs=%d, number of coupled chains per run=%d => total of %d chains \n", gAInfo.numberOfRuns, gAInfo.numberCoupledChains, totalNumChains ); 

#ifdef DEBUG_LNL_VERIFY
  gAInfo.debugTree = new TreeAln( byteFileName); 
#endif

  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    {
      // initialize some trees 
      vector<TreeAln*> trees; 
      if(i == 0 )	
	{
	  trees.push_back(new TreeAln(byteFileName)); 
	}
      else 
	{
	  CoupledChains* aRun = runs[0]; 
	  for(int i = 0; i < aRun->getNumberOfChains(); ++i)	    
	    {
	      TreeAln *traln = aRun->getChain(i)->traln; 
	      trees.push_back(new TreeAln(*traln)); 
	    }	  
	}      
    }

  
  // only use one restorer for all chains 
  LnlRestorer *restorer = new LnlRestorer(runs[0]->getChain(0));
  for(auto run : runs )
    {
      for(int i = 0; i < run->getNumberOfChains(); ++i )
	{
	  Chain *chain = run->getChain(i); 
	  chain->setRestorer(restorer); 
	}
    }



  // initialize following runs 

  assert(0); 
  // TODO all of this must be integrated 
#if 0 
  

  for(int i = 0; i < totalNumChains; ++i)
    {       
      tree *myTree = theChain->traln->getTr();       

      setupProposals(theChain, initParams); 

      /* init the param dump  */      
      initParamDump(traln, &(theChain->dump)); 

      // for(int j = 0; j < traln->getNumberOfPartitions(); ++j ) 
      // 	theChain->traln->initRevMat(j);

      // initLocalRng(theChain); 
      
      if( i % gAInfo.numberCoupledChains == 0)
	{
	  /* initialize with new tree  */
	  if( i  / gAInfo.numberCoupledChains < gAInfo.numberOfStartingTrees )
	    {
	      boolean hasBranchLength = readTreeWithOrWithoutBL(theChain->traln->getTr(), treeFH); 

	      if(hasBranchLength)
		{
		  int count = 0; 
		  tree *tr = theChain->traln->getTr();
		  traverseInitCorrect(tr->start->back, &count, theChain->traln ) ;
		  assert(count == 2 * theChain->traln->getTr()->mxtips -3); 
		}
	      else 
		{
		  /* set some standard branch lengths */
		  /* TODO based on prior?   */
		  int count = 0; 
		  traverseInitFixedBL( theChain->traln->getTr()->start->back, &count, theChain->traln, INIT_BRANCH_LENGTHS );
		  assert(count  == 2 * theChain->traln->getTr()->mxtips - 3);	  	      
		}

	      PRINT("initializing chain %d with provided starting tree\n", i); 
	    }
	  else 
	    {
	      exa_makeRandomTree(theChain->traln->getTr());
	      PRINT("initializing chain %d with random tree\n", i); 

	      /* TODO maybe prior for initial branch lengths   */
	      int count = 0; 
	      traverseInitFixedBL( theChain->traln->getTr()->start->back, &count, theChain->traln, INIT_BRANCH_LENGTHS );

	      assert(count  == 2 * theChain->traln->getTr()->mxtips - 3);	  
	    }
	}
      else 
	{
	  Chain *coldChain = *runs + (theChain->id / gAInfo.numberCoupledChains) * gAInfo.numberCoupledChains; 
	  copyState(theChain , coldChain);
	  applyChainStateToTree(theChain); 
	}

      // evaluateFullNoBackup(theChain);

      // /* now save the tree to a chain chains */
      // saveTreeStateToChain(theChain); 

      // PRINT("Initial LnL for chain %d is  %f\ttree-length=%.3f\tseed=%u,%u\n", theChain->id, theChain->traln->getTr()->likelihood, 
      // 	    branchLengthToReal(theChain->traln->getTr(), getTreeLength(theChain->traln, theChain->traln->getTr()->nodep[1]->back)),
      // 	    theChain->rKey.v[0], theChain->rKey.v[1]
      // 	    ); 

      if(isOutputProcess() )
	{	  
	  if( i % gAInfo.numberCoupledChains == 0)
	    {
	      makeChainFileNames(theChain, i / gAInfo.numberCoupledChains); 
	      // printf("\n\ninitializing output file\n\n") ;
	      initializeOutputFiles(theChain); 
	    }
	  else 
	    {	      
	      Chain *previousChain = *runs + (i  - 1 ) ; 
	      theChain->topologyFile = previousChain->topologyFile; 
	      theChain->outputParamFile = previousChain->topologyFile; 
	    }
	}
    }


#endif 


  {
    Chain* chain = runs[0]->getChain(0); 
    int numTax = chain->traln->getTr()->mxtips; 
    gAInfo.bipHash = new BipartitionHash(numTax, gAInfo.numberOfRuns);
  }
  
  
  if(gAInfo.numberOfStartingTrees > 0)
    fclose(treeFH); 


  PRINT("\n"); 

  // TODO kill initial tree 
}



void traverseInitFixedBL(nodeptr p, int *count, TreeAln *traln,  double z )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;
  
  for( i = 0; i < traln->getNumBranches(); i++)
      traln->setBranchLengthSave(z, i, p); 
  
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitFixedBL(q->back, count, traln, z);
	  q = q->next;
	} 
    }
}


void traverseAndTreatBL(node *p, TreeAln *traln, double *blBuf, int* cnt, boolean restore)
{
  tree *tr = traln->getTr();
  nodeptr q; 
  assert(traln->getNumBranches() == 1); 

  if(restore == TOPO_RESTORE )
    traln->setBranchLengthSave(blBuf[*cnt], 0,p); 
  else if(restore == TOPO_SAVE)    
    blBuf[*cnt] =  traln->getBranchLength(p,0); 
  else 
    assert(0); 
  *cnt += 1 ; 
  
  if( NOT isTip(p->number, tr->mxtips))
    {
      q = p->next; 
      while(q != p )
	{
	  traverseAndTreatBL(q->back, traln, blBuf, cnt, restore); 
	  q = q->next; 
	}
    }
}





/**
   @brief Applies the Chain of the chain to its tree. 

   Notice: you must not simply change the tree pointer. More
   modifications are necessary to do so.

   @param boolean checkLnl -- should we check, if the lnl is the same as
   before? If we applied it the first time, there is no before.
 */ 
void applyChainStateToTree(Chain *chain)
{
  tree *tr = chain->traln->getTr(); 
  
  /* TODO enable multi-branch    */
  assert(chain->traln->getNumBranches() == 1); 

  chain->dump.topology->restoreTopology(*(chain->traln));

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, chain->traln, chain->dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 

  /* restore model parameters */
  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i ; 

      pInfo *partition = chain->traln->getPartition( i);
      partition->alpha = info->alpha;

      memcpy(partition->substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(partition->frequencies, info->frequencies, 4 * sizeof(double));
      
      chain->traln->initRevMat(i);
      chain->traln->discretizeGamma(i);	 
    } 

  evaluateFullNoBackup(chain); 

  chain->wasAccepted = TRUE; 

}



/**
   @brief Save all relevan information from the tree into the chain
   Chain. 
 */ 
void saveTreeStateToChain(Chain *chain)
{
  tree *tr  = chain->traln->getTr();     
  chain->dump.topology->saveTopology(*(chain->traln));

  /* save branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, chain->traln, chain->dump.branchLengths, &cnt, TOPO_SAVE); 
  assert(cnt == 2 * tr->mxtips - 3 ); 

  /* save model parameters */
  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i; 
      pInfo *partition = chain->traln->getPartition(i); 
      
      info->alpha =  partition->alpha; 
      
      memcpy(info->substRates, partition->substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, partition->frequencies, 4 * sizeof(double)); 
    }  


}





/**
   @brief draws a proposal function.

   Notice: this could be extended later, if we decide to make this
   dependent on the previous state.
   
   Furthermore, we must be sure now that category weights and relative
   proposal weights sum up to 1 each. 
   
 */ 
void drawProposalFunction(Chain *chain, proposalFunction **result )
{  
  
  *result = NULL; 
  category_t
    cat = category_t(drawSampleProportionally(chain,chain->categoryWeights, NUM_PROP_CATS) + 1); /* it is 1-based */

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

