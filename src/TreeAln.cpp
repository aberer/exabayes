
#include <sstream>
#include <cassert>

#include "config.h"
#include "TreeAln.hpp"

#include "GlobalVariables.hpp"
#include "output.h"
#include "branch.h"


// here we initialize the max/min values for our various
// parameters. Static const is essentially like a global variable but
// saver. You access it with e.g., TreeAln::zmin, since the variable
// does not belong to an instance of the class, but the class in general. 
// It is a bit over-engineering. 

// most of it has been copied over from raxml. But maybe we want to differ? 



// TODO these raxml values would allow for branch lengths up to 13,
// I'll reduce it now because i ran into trouble there, but this
// definitely needs to be discussed
// #if 0 
const double TreeAln::zMin = 1.0E-15 ; 
const double TreeAln::zMax = (1.0 - 1.0E-6) ; 
// #else 
// const double TreeAln::zMin = 1.0E-8 ; 
// const double TreeAln::zMax = (1.0 - 1.0E-3) ; 
// #endif


const double TreeAln::zZero = TreeAln::zMax + ( 1 - TreeAln::zMax) / 2 ; 

const double TreeAln::rateMin = 0.0000001; 
const double TreeAln::rateMax = 1000000.0; 

const double TreeAln::alphaMin = 0.02; 
const double TreeAln::alphaMax = 1000.0; 

const double TreeAln::freqMin = 0.001; 


const double TreeAln::initBL = 0.65; // TODO I'd prefer absolute real value of 0.1  (currentyl 0.15)


TreeAln::TreeAln()
  : parsimonyEnabled(true)
{
  tree *tre = (tree*)exa_calloc(1,sizeof(tree));
  this->tr = tre; 
  initDefault();
}



TreeAln::~TreeAln()
{
  if(tr->aliaswgt != NULL)
    exa_free(tr->aliaswgt); 
  if(tr->rateCategory != NULL)
    exa_free(tr->rateCategory); 
  if(tr->patrat != NULL)
    exa_free(tr->patrat);
  if(tr->patratStored != NULL)
    exa_free(tr->patratStored);
  if(tr->lhs != NULL)
    exa_free(tr->lhs);  
  if(tr->yVector != NULL)
    exa_free(tr->yVector); 
  if(tr->nameList != NULL)
    {
      for(int i = 1; i <= tr->mxtips; ++i )
	exa_free(tr->nameList[i]);
      exa_free(tr->nameList);
    }

  if(tr->tree_string != NULL)
    exa_free(tr->tree_string);
  if(tr->tree0 != NULL)
    exa_free(tr->tree0);
  if(tr->tree1 != NULL)
    exa_free(tr->tree1);
  
  if(tr->constraintVector != NULL)
    exa_free(tr->constraintVector); 
  if(tr->nodeBaseAddress != NULL)
    exa_free(tr->nodeBaseAddress);

  
  int numPart = getNumberOfPartitions();
  for(int i = 0; i < numPart ;++i)
    {
      pInfo *partition = getPartition(i); 
      exa_free(partition); 
    }
  
  if(numPart > 0)
    {
#if HAVE_PLL != 0 
      if(getNumberOfPartitions() > 0 )
	exa_free(partitions->partitionData); 
      exa_free(partitions); 
#else 
      exa_free(tr->partitionData);
#endif
    }

  if(parsimonyEnabled)
#if HAVE_PLL != 0
    freeParsimonyDataStructures(tr, partitions);  
#else 
  ; 
#endif

  exa_free(tr);
}


void TreeAln::initializeFromByteFile(string _byteFileName)
{
#if HAVE_PLL != 0
  partitionList *pl = (partitionList*)exa_calloc(1,sizeof(partitionList)); 
  partitions = pl;
  this->initializeTreePLL(_byteFileName);

#else 
  extern char byteFileName[1024]; 
  strcpy(byteFileName, _byteFileName.c_str() ) ;  

  analdef adef ; 

  adef.max_rearrange          = 21;
  adef.stepwidth              = 5;
  adef.initial                = 10;
  adef.bestTrav               = 10;
  adef.initialSet             = FALSE; 
  adef.mode                   = BIG_RAPID_MODE; 
  adef.likelihoodEpsilon      = 0.1; 
  adef.permuteTreeoptimize    = FALSE; 
  adef.perGeneBranchLengths   = FALSE;   
  adef.useCheckpoint          = FALSE;
  
  initializeTree(tr, &adef);   
#endif  

  for(int i = 0; i < getNumberOfPartitions(); ++i)
    {

      // TODO empirical freqs? 

      // TODO aa? 

      pInfo *partition =  getPartition(i);

      assert(partition->dataType == DNA_DATA);

      vector<double> tmp(partition->states, 1.0 / (double)partition->states ); 
      setFrequenciesBounded(tmp,i);
      vector<double>tmp2( numStateToNumInTriangleMatrix(partition->states), 1.0 / (double) numStateToNumInTriangleMatrix(partition->states)); 
      setRevMatBounded(tmp2,i);
    }
}

double TreeAln::getTreeLengthExpensive() const
{
  vector<branch> branches; 
  extractBranches(*this, branches); 
  
  assert(getNumBranches() == 1 ); 

  double result = 1; 
  for(auto b : branches)
    result *= b.length[0]; 

  return result; 
}

void TreeAln::verifyTreeLength() const
{
#if TODO 
  double tlVerified = getTreeLengthExpensive() ;
  assert(treeLength == tlVerified); 
#endif
}

void TreeAln::enableParsimony()
{
#if HAVE_PLL == 0
  // allocateParsimonyDataStructures(tr);   
  // cout << "parsimony not implemented yet in ExaML" << endl; 
  // assert(0); 
#else 
  allocateParsimonyDataStructures(tr, partitions);   
#endif
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




void TreeAln::clipNodeDefault(nodeptr p, nodeptr q )
{
  double tmp = TreeAln::initBL; 
  clipNode(p,q, tmp);
}


void TreeAln::clipNode(nodeptr p, nodeptr q, double &z)
{
  p->back = q; 
  q->back = p; 
  assert(getNumBranches() == 1);   

  setBranchLengthBounded(z,0,p);
  setBranchLengthBounded(z,0,q);
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

#define UNSAFE_EXACT_TREE_COPY

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
      memcpy(partitionLhs->frequencies, partitionRhs->frequencies, partitionRhs->states * sizeof(double)); 
      memcpy(partitionLhs->substRates, partitionRhs->substRates, numStateToNumInTriangleMatrix(partitionRhs->states) * sizeof(double)); 
      partitionLhs->alpha = partitionRhs->alpha; 
      this->initRevMat(i);
      this->discretizeGamma(i);       
    }

  
  this->unlinkTree();
  int mxtips = rhs.getTr()->mxtips ;
  tree *rhsTree = rhs.getTr(),
    *thisTree = getTr();


#ifdef UNSAFE_EXACT_TREE_COPY
  // if this works, it is one of the most hackney things, I've ever done... 
  for(int i = mxtips+1 ; i < 2* tr->mxtips-1; ++i)
    {
      nodeptr rhsNode = rhsTree->nodep[i],
	lhsNode = thisTree->nodep[i];       
      hookup(lhsNode, lhsNode -   (rhsNode - rhsNode->back), rhsNode->z, getNumBranches()) ;             

      lhsNode = lhsNode->next; 
      rhsNode = rhsNode->next; 
      hookup(lhsNode, lhsNode -   (rhsNode - rhsNode->back), rhsNode->z, getNumBranches()) ;       

      lhsNode = lhsNode->next; 
      rhsNode = rhsNode->next; 
      hookup(lhsNode, lhsNode  -   (rhsNode - rhsNode->back), rhsNode->z, getNumBranches()) ; 
    }
#else 


  for(int i = 1 ; i <= mxtips; ++i )
    {
      nodeptr a = rhsTree->nodep[i],
	b = rhsTree->nodep[i]->back; 

#ifdef DEBUG_LINK_INFO
      cout << "hooking up tip " << a->number << " and " << b->number << endl; 
#endif
      
      hookup(this->getUnhookedNode(a->number), this->getUnhookedNode(b->number),a->z, getNumBranches()); 
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

#endif

  tr->start = tr->nodep[rhsTree->start->number]; 
  debug_checkTreeConsistency(*this);
  
  return *this; 
}



int TreeAln::getNumBranches() const
{
#if HAVE_PLL != 0
  return partitions->perGeneBranchLengths  ? partitions->numberOfPartitions : 1 ; 
#else 
  return tr->numBranches; 
#endif
}



int TreeAln::getNumberOfPartitions() const
{
#if HAVE_PLL != 0
  return partitions->numberOfPartitions; 
#else 
  return tr->NumberOfModels; 
#endif
}



// usefull stuff 
pInfo* TreeAln::getPartition(int model)  const
{
  assert(model < getNumberOfPartitions()); 
  
#if HAVE_PLL != 0  
  return partitions->partitionData[model]; 
#else 
  return tr->partitionData + model; 
#endif
}





// this is BS 
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



void TreeAln::setFrequenciesBounded(vector<double> &newValues, int model )
{  
  // cout << "set freq:"; 
  // for(auto &v : newValues)
  //   cout << v << "," ; 
  // cout << endl; 

  int changed = 0; 
  double sum =0 ; 
  for(double &v : newValues)
    if(v < freqMin)
      {
	v = freqMin; 
	changed++; 
      }
    else 
      sum += v; 

  sum += ( changed * freqMin ); 
  for(double &v : newValues)
    if(v != freqMin)
      v /= sum; 


  // cout << "after norm:"; 
  // for(auto &v : newValues)
  //   cout << v << "," ; 
  // cout << endl; 	  

  pInfo *partition = getPartition(model); 
  for(nat i = 0; i < newValues.size(); ++i)
    partition->frequencies[i] = newValues[i];   
  
  sum = 0; 
  for(auto &v : newValues)
    sum += v; 
  assert(fabs (sum - 1.0) < 1e-6  );  
}



void TreeAln::setRevMatBounded(vector<double> &newValues, int model)
{ 
  // cout << "setting revmat "; 
  // for(auto &v : newValues)
  //   cout << v << ","; 
  // cout << endl; 

  vector<double> normValues = newValues; 

  for(auto &v : normValues)
    v /= normValues[normValues.size()-1]; 
  
  double sum = 0; 
  for(double &v : normValues)
    {
      if(v < rateMin)
	v = rateMin; 
      else if(rateMax < v )
	v = rateMax; 
      sum += v ; 
    }

  // cout << "after correction " ; 
  // for(auto &v : normValues)
  //   cout << v << "," ; 
  // cout << endl; 

  pInfo *partition = getPartition(model); 
  for(nat i = 0; i < newValues.size(); ++i)
    partition->substRates[i] = newValues[i]; 
}


/** 
    @brief save setting method for a branch length

    @param newValue -- INTERNAL branch length
 */  
void TreeAln::setBranchLengthBounded(double &newValue, int model, nodeptr p)
{
#if TODO 
  double oldZ = p->z[model] ; 
#endif
  if(newValue < zMin)
    newValue = zMin; 
  if (zMax < newValue)
    newValue = zMax; 

#if TODO 
  if(newValue != oldZ)
    treeLength *= newValue / oldZ ; 
#endif
  
  p->z[model] = p->back->z[model] = newValue; 
}


/** 
    @brief save setting method for an alpha value 
    @return newValue after check 
    
 */  
void TreeAln::setAlphaBounded(double &newValue, int model)
{
  if(newValue < alphaMin)
    newValue = alphaMin; 
  if(alphaMax < newValue )
    newValue = alphaMax; 
  
  pInfo *partition = getPartition(model);
  partition->alpha =  newValue; 
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

void TreeAln::initializeTreePLL(string byteFileName)
{
  partitionList *pl = getPartitionsPtr();
  pl->partitionData = (pInfo**)exa_malloc(NUM_BRANCHES*sizeof(pInfo*));

  double **empFreq = NULL; 
  bool perGeneBL = false; 
  initializePartitionsPLL(byteFileName, &empFreq, perGeneBL );

  initializePartitionsSequential(getTr(), getPartitionsPtr());
  initModel(getTr(), empFreq, getPartitionsPtr());
} 


void TreeAln::initializePartitionsPLL(string byteFileName, double ***empiricalFrequencies, bool multiBranch)
{
  tree *tr = getTr();
  partitionList *partitions = getPartitionsPtr();

  unsigned char *y;

  FILE 
    *byteFile = myfopen(byteFileName.c_str(), "rb");	 

  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
  myBinFread(&(partitions->numberOfPartitions),  sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);


  partitions->perGeneBranchLengths = multiBranch ? 1 : 0 ;

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


ostream& operator<< (ostream& out,  TreeAln&  traln)
{
  TreePrinter tp(true, false, false); 
  return out << tp.printTree(traln); 
}



void TreeAln::collapseBranch(branch b)
{
  assert(getNumBranches() == 1 ); 
  nodeptr p = findNodeFromBranch(  getTr(), b); 
  p->z[0] = p->back->z[0] = TreeAln::zZero;   
}


bool TreeAln::isCollapsed(branch b ) 
{
  assert(getNumBranches() == 1 ); 
  nodeptr p = findNodeFromBranch(  getTr(), b); 
  return p->z[0] >=  TreeAln::zMax ; 
}


void TreeAln::setBranchLengthUnsafe(branch b ) 
{
  assert(getNumBranches() == 1 ); 
  nodeptr p = findNodeFromBranch(  getTr(), b);   
  p->z[0] = p->back->z[0] = b.length[0];   
} 


vector<double> TreeAln::getRevMat(int model) const 
{
  vector<double> result; 
  pInfo *partition = getPartition(model); 
  double sum = 0; 
  for(int i =0 ; i < numStateToNumInTriangleMatrix(partition->states); ++i)
    {
      result.push_back(partition->substRates[i]);
      sum += result[i]; 
    }

  for_each(result.begin(), result.end(), [&](double &v) { v /= sum ; }) ; 

  return result; 
}

vector<double> TreeAln::getFrequencies(int model) const
{
  vector<double> result; 
  pInfo* partition = getPartition(model) ; 
  for(int i = 0; i < partition->states; ++i) 
    result.push_back(partition->frequencies[i]); 
  return result; 
}

bool TreeAln::revMatIsImmutable(int model) const
{
#ifdef UNSURE
  assert(0); 
#endif
    
  pInfo *partition = getPartition(model); 
    
  return partition->states == 20 && partition->protModels != GTR; 
} 











// bool TreeAln::toString(bool showInternal, bool showBL, bool printNames) const
// {

//   assert(not printNames); 

//   tree *tr = traln.getTr();
//   Tree2stringNexus(tr->tree_string, tr, tr->start->back,0);
//   stringstream ss; 
//   ss << tr->tree_string;   
//   return ss.str(); 
//   return out; 
// }


