
#include "config.h"
#include "TreeAln.hpp"

#include "globals.h"
#include "main-common.h" 
#include "output.h"



// #define DEBUG_LINK_INFO 	// TODO erase 


/** 
    @param tr -- tree to be initialized from 
 */ 
TreeAln::TreeAln(tree *trOuter) 
{  
  tree *newTr = (tree*)exa_calloc(1, sizeof(tree) ); 
  newTr->mxtips = trOuter->mxtips; 
  this->tr = newTr; 
  initDefault(); 

  // create a new tree 
#if  ( HAVE_PLL == 1  ) 
  partitionList *partition = (partitionList*)exa_calloc(1,sizeof(partitionList)); 
  this->partitions = partition; 
  initializeTree(newTr, partition ,gAInfo.adef); 
  this->partitions = partition; 
#else 
  initializeTree(newTr ,gAInfo.adef); 
#endif        
  initFromTree(trOuter);  
}


void TreeAln::initFromTree(tree *trOuter) 
{
  for(int j = 1; j <= trOuter->mxtips; ++j)
    this->tr->nodep[j]->hash = trOuter->nodep[j]->hash; 
  this->tr->bitVectors = trOuter->bitVectors; 
}



/** 
    @brief standard values for the tree 
*/ 
void TreeAln::initDefault()
{   
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA; 
  tr->multiStateModel  = GTR_MULTI_STATE;
#if (HAVE_PLL == 0 ) 
  tr->useGappedImplementation = FALSE;
  tr->saveBestTrees          = 0;
#endif
  tr->saveMemory = FALSE;
  tr->manyPartitions = FALSE;
  tr->categories             = 25;
  tr->grouped = FALSE;
  tr->constrained = FALSE;
  tr->gapyness               = 0.0; 
  tr->useMedian = FALSE;
}



TreeAln::~TreeAln()
{
  // TODO 
}


void TreeAln::unlinkTree()
{
#ifdef DEBUG_LINK_INFO
  cout << "unlinking everything" << endl; 
#endif

  // unlink tips 
  for(int i = 1; i < tr->mxtips+1; ++i)
    {
      nodeptr p = tr->nodep[i]; 
      p->back = NULL; 
    }

  for(int i = tr->mxtips + 1; i < 2 * tr->mxtips; ++i)
    {      
      nodeptr p = tr->nodep[i]; 
      p->back = NULL; 
      p->next->back = NULL; 
      p->next->next->back = NULL; 
    }

}


nodeptr TreeAln::getUnhookedNode(int number)
{  
  tree *tr = getTr();
  nodeptr p = tr->nodep[number]; 

  if(isTip(number, tr->mxtips) )
    {
      assert(p->back == NULL); 
      return p; 
    }


  nodeptr q = p ; 
  do 
    {
      if(q->back == NULL)
	return q ; 
      q = q->next; 
    } while(p != q); 

  assert(0);
  return NULL;
}



/**
   @brief copies the entire state from the rhs to this tree/alignment.
   
   all parameters are copied and initialized, topology and branch
   lengths are copied.
 */ 
TreeAln& TreeAln::operator=( TreeAln& rhs)
{
  assert(&rhs != this); 

  // copy partition parameters 
  for(int i = 0; i < rhs.getNumberOfPartitions(); ++i)
    {
      pInfo *partitionRhs = rhs.getPartition(i); 
      pInfo *partitionLhs = this->getPartition(i); 
      memcpy(partitionLhs->frequencies, partitionRhs->frequencies, 4 * sizeof(double)); 
      memcpy(partitionLhs->substRates, partitionRhs->substRates, 6 * sizeof(double)); 
      partitionLhs->alpha = partitionRhs->alpha; 
      this->initRevMat(i);
      makeGammaCats(partitionLhs->alpha, partitionLhs->gammaRates, 4, tr->useMedian);
    }
  

  this->unlinkTree();
  int mxtips = rhs.getTr()->mxtips ;
  tree *rhsTree = rhs.getTr(),
    *thisTree = getTr();
  for(int i = 1 ; i <= mxtips; ++i )
    {
      nodeptr a = rhsTree->nodep[i],
	b = rhsTree->nodep[i]->back; 

#ifdef DEBUG_LINK_INFO
      cout << "hooking up tip " << a->number << " and " << b->number << endl; 
#endif
      
      hookup(this->getUnhookedNode(a->number), 
	     this->getUnhookedNode(b->number),a->z, getNumBranches()); 
    }

  for(int i = mxtips+1; i < 2 * tr->mxtips-1; ++i)
    {      
      nodeptr pRhs = rhsTree->nodep[i], 
	q = pRhs; 
      do 
	{
	  if(NOT branchExists(thisTree, constructBranch(q->number, q->back->number)))
	    {
#ifdef DEBUG_LINK_INFO
	      cout << "hooing up " << q->number << " and " << q->back->number << endl; 
#endif
	      hookup(getUnhookedNode(q->number ) , getUnhookedNode(q->back->number) , q->z, getNumBranches()); 
	    }
	  else 
	    {
#ifdef DEBUG_LINK_INFO
	      cout << "branch between " << q->number <<  " and "  << q->back->number << " already existing" <<endl; 
#endif
	    }

	  q = q->next ; 
	} while(q != pRhs); 
    }

  tr->start = tr->nodep[rhsTree->start->number]; 
  debug_checkTreeConsistency(this->tr);
  
  return *this; 
}



/**
   @brief shallow copy of the tree 
 */ 
TreeAln::TreeAln(const TreeAln& rhs)  
{  
  this->tr = rhs.tr; 
#if HAVE_PLL == 1 
  this->partitions = rhs.partitions; 
#endif
}



int TreeAln::getNumBranches()
{
#if HAVE_PLL == 1 
  return partitions->perGeneBranchLengths  ? partitions->numberOfPartitions : 1 ; 
#else 
  return tr->numBranches; 
#endif
}



int TreeAln::getNumberOfPartitions()
{
#if HAVE_PLL == 1 
  return partitions->numberOfPartitions; 
#else 
  return tr->NumberOfModels; 
#endif
}



// usefull stuff 
pInfo* TreeAln::getPartition(int model) 
{
#if HAVE_PLL == 1   
  return partitions->partitionData[model]; 
#else 
  return tr->partitionData + model; 
#endif
}


int& TreeAln::accessExecModel(int model)
{
#if HAVE_PLL == 1 
  return partitions->partitionData[model]->executeModel; 
#else 
  return tr->executeModel[model]; 
#endif
}


double& TreeAln::accessPartitionLH(int model)
{
#if HAVE_PLL == 1 
  return partitions->partitionData[model]->partitionLH; 
#else  
  return tr->perPartitionLH[model]; 
#endif
}


void TreeAln::initRevMat(int model)
{
#if HAVE_PLL == 1
  initReversibleGTR(tr, getPartitionsPtr() , model); 
#else 
  initReversibleGTR(tr, model); 
#endif
}




