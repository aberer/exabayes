

#include "config.h"
#include "TreeAln.hpp"

#include "globals.h"
#include "main-common.h" 
#include "output.h"


// here we initialize the max/min values for our various
// parameters. Static const is essentially like a global variable but
// saver. You access it with e.g., TreeAln::zmin, since the variable
// does not belong to an instance of the class, but the class in general. 
// It is a bit over-engineering. 

// most of it has been copied over from raxml. But maybe we want to differ? 

const double TreeAln::zMin = 1.0E-15 ; 
const double TreeAln::zMax = (1.0 - 1.0E-6) ; 

const double TreeAln::rateMin = 0.0000001; 
const double TreeAln::rateMax = 1000000.0; 

const double TreeAln::alphaMin = 0.02; 
const double TreeAln::alphaMax = 1000.0; 

const double TreeAln::freqMin = 0.001; 


// #define DEBUG_LINK_INFO 	// TODO erase 


TreeAln::TreeAln(char *bytefile)
{ 
  tree *tre = (tree*)exa_calloc(1,sizeof(tree)) ; 
  this->tr = shared_ptr<tree>(tre); 

#if HAVE_PLL != 0
  partitionList *pl = (partitionList*)exa_calloc(1,sizeof(partitionList)); 
  partitions = shared_ptr<partitionList>(pl);
  this->initializeTreePLL();
#else 
  initializeTree(tre, gAInfo.adef);   
#endif
}

TreeAln::TreeAln(const TreeAln &rhs)
  : tr(rhs.tr)
#if HAVE_PLL != 0
  , partitions(rhs.partitions)
#endif
{  
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
#if HAVE_PLL == 0 
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
      this->discretizeGamma(i);       
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
  debug_checkTreeConsistency(this->getTr());
  
  return *this; 
}



// /**
//    @brief shallow copy of the tree 
//  */ 
// TreeAln::TreeAln(const TreeAln& rhs)  
// {  
//   this->tr = rhs.tr; 
// #if HAVE_PLL != 0
//   this->partitions = rhs.partitions; 
// #endif
// }



int TreeAln::getNumBranches()
{
#if HAVE_PLL != 0
  return partitions->perGeneBranchLengths  ? partitions->numberOfPartitions : 1 ; 
#else 
  return tr->numBranches; 
#endif
}



int TreeAln::getNumberOfPartitions()
{
#if HAVE_PLL != 0
  return partitions->numberOfPartitions; 
#else 
  return tr->NumberOfModels; 
#endif
}



// usefull stuff 
pInfo* TreeAln::getPartition(int model) 
{
  assert(model < getNumberOfPartitions()); 
  
#if HAVE_PLL != 0  
  return partitions->partitionData[model]; 
#else 
  return tr->partitionData + model; 
#endif
}


int& TreeAln::accessExecModel(int model)
{
#if HAVE_PLL != 0
  return partitions->partitionData[model]->executeModel; 
#else 
  return tr->executeModel[model]; 
#endif
}


double& TreeAln::accessPartitionLH(int model)
{
#if HAVE_PLL != 0
  return partitions->partitionData[model]->partitionLH; 
#else  
  return tr->perPartitionLH[model]; 
#endif
}


void TreeAln::initRevMat(int model)
{
#if HAVE_PLL != 0
  initReversibleGTR(getTr(), getPartitionsPtr() , model); 
#else 
  initReversibleGTR(getTr(), model); 
#endif
}




/** 
    @brief save setting method for a frequency parameter 

    @return newValue after check
 */  
double TreeAln::setFrequencySave(double newValue, int model, int position )
{
  if(newValue < freqMin )
    newValue = freqMin; 

  pInfo *partition = getPartition(model); 
  // TODO assert? 
  partition->frequencies[position] = newValue; 
  return newValue; 
}

/** 
    @brief save setting method for a substitution parameter 
    @return newValue after check
 */  
double TreeAln::setSubstSave(double newValue, int model, int position)
{
  if(newValue < rateMin)
    newValue = rateMin; 
  if(rateMax < newValue )
    newValue = rateMax; 

  pInfo *partition = getPartition(model);
  assert(position < partition->states * partition->states - partition->states); 
  partition->substRates[position] = newValue; 
  
  return newValue; 
}



/** 
    @brief save setting method for a branch length
 */  
double TreeAln::setBranchLengthSave(double newValue, int model, nodeptr p)
{
  if(newValue < zMin)
    newValue = zMin; 
  if (zMax < newValue)
    newValue = zMax; 
  
  p->z[model] = p->back->z[model] = newValue; 

  return newValue; 
}


/** 
    @brief save setting method for an alpha value 
    @return newValue after check 
    
 */  
double TreeAln::setAlphaSave(double newValue, int model)
{
  if(newValue < alphaMin)
    newValue = alphaMin; 
  if(alphaMax < newValue )
    newValue = alphaMax; 

  pInfo *partition = getPartition(model);
  partition->alpha =  newValue; 
  
  return newValue; 
}


/**
   @brief makes the discrete categories for the gamma
   distribution. Has to be called, if alpha was updated.   

 */ 
void TreeAln::discretizeGamma(int model)
{
  pInfo *partition =  getPartition(model); 
  makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);
}




// initialization code from the PLL in future version, we could also
// support plain files again. Just do not really know what to do with
// it =/
#if HAVE_PLL  != 0 

void TreeAln::initializeTreePLL()
{
  partitionList *pl = getPartitionsPtr();
  pl->partitionData = (pInfo**)exa_malloc(NUM_BRANCHES*sizeof(pInfo*));

  double **empFreq = NULL; 
  initializePartitionsPLL(&(byteFileName[0]), &empFreq);

  initializePartitionsSequential(getTr(), getPartitionsPtr());
  initModel(getTr(), empFreq, getPartitionsPtr());
} 


void TreeAln::initializePartitionsPLL(char *bytefile, double ***empiricalFrequencies)
{
  analdef *adef = gAInfo.adef ; 
  tree *tr = getTr();
  partitionList *partitions = getPartitionsPtr();

  unsigned char *y;

  FILE 
    *byteFile = myfopen(bytefile, "rb");	 

  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
  myBinFread(&(partitions->numberOfPartitions),  sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);


  partitions->perGeneBranchLengths = adef->perGeneBranchLengths;

  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed RAxML-specific */

  tr->aliaswgt                   = (int *)exa_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
  myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       

  tr->rateCategory    = (int *)    exa_malloc((size_t)tr->originalCrunchedLength * sizeof(int));	  

  tr->patrat          = (double*)  exa_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->patratStored    = (double*)  exa_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
  tr->lhs             = (double*)  exa_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 

  *empiricalFrequencies = (double **)exa_malloc(sizeof(double *) * partitions->numberOfPartitions);

  y = (unsigned char *)exa_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));
  tr->yVector = (unsigned char **)exa_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

  for( nat i = 1; i <= (size_t)tr->mxtips; i++)
    tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength]; 

  setupTree(tr, FALSE, partitions);

  for(int i = 0; i < partitions->numberOfPartitions; i++)
    partitions->partitionData[i]->executeModel = TRUE;

  /* data structures for convergence criterion need to be initialized after! setupTree */

  for( nat i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int len;
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr->nameList[i] = (char*)exa_malloc(sizeof(char) * (size_t)len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      /*printf("%s \n", tr->nameList[i]);*/
    }  

  for( nat i = 1; i <= (size_t)tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  for(int model = 0; model < partitions->numberOfPartitions; model++)
    {
      int 
	len;

      pInfo 
	*p = partitions->partitionData[model];

      myBinFread(&(p->states),             sizeof(int), 1, byteFile);
      myBinFread(&(p->maxTipStates),       sizeof(int), 1, byteFile);
      myBinFread(&(p->lower),              sizeof(int), 1, byteFile);
      myBinFread(&(p->upper),              sizeof(int), 1, byteFile);
      myBinFread(&(p->width),              sizeof(int), 1, byteFile);
      myBinFread(&(p->dataType),           sizeof(int), 1, byteFile);
      myBinFread(&(p->protModels),         sizeof(int), 1, byteFile);
      myBinFread(&(p->autoProtModels),     sizeof(int), 1, byteFile);
      myBinFread(&(p->protFreqs),          sizeof(int), 1, byteFile);
      myBinFread(&(p->nonGTR),             sizeof(boolean), 1, byteFile);
      myBinFread(&(p->numberOfCategories), sizeof(int), 1, byteFile); 

      /* later on if adding secondary structure data

	 int    *symmetryVector;
	 int    *frequencyGrouping;
      */

      myBinFread(&len, sizeof(int), 1, byteFile);
      p->partitionName = (char*)exa_malloc(sizeof(char) * (size_t)len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);

      (*empiricalFrequencies)[model] = (double *)exa_malloc(sizeof(double) * (size_t)partitions->partitionData[model]->states);
      myBinFread((*empiricalFrequencies)[model], sizeof(double), partitions->partitionData[model]->states, byteFile);
    }

  myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

  fclose(byteFile);
}


#endif
