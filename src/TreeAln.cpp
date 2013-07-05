#include <sstream>
#include <cassert>
#include <cstring>

#include "axml.h"

#include "TreeAln.hpp"
#include "GlobalVariables.hpp"

#include "BoundsChecker.hpp"

const double TreeAln::zZero = BoundsChecker::zMax + ( 1 - BoundsChecker::zMax) / 2 ; 
const double TreeAln::initBL = 0.65; // TODO I'd prefer absolute real value of 0.1  (currentyl 0.15)

TreeAln::TreeAln()
  : parsimonyEnabled(true)
{
  memset(&tr,0,sizeof(tree)); 
  initDefault();
}


TreeAln::~TreeAln()
{  
  // this is slightly crazy, but we'll have to do the cleanup
  // somewhere.

  if(tr.aliaswgt != NULL)
    exa_free(tr.aliaswgt); 
  if(tr.rateCategory != NULL)
    exa_free(tr.rateCategory); 
  if(tr.patrat != NULL)
    exa_free(tr.patrat);
  if(tr.patratStored != NULL)
    exa_free(tr.patratStored);
  if(tr.lhs != NULL)
    exa_free(tr.lhs);  
  if(tr.yVector != NULL)
    exa_free(tr.yVector); 
  if(tr.nameList != NULL)
    {
      for(int i = 1; i <= tr.mxtips; ++i )
	exa_free(tr.nameList[i]);
      exa_free(tr.nameList);
    }

  if(tr.tree_string != NULL)
    exa_free(tr.tree_string);
  if(tr.tree0 != NULL)
    exa_free(tr.tree0);
  if(tr.tree1 != NULL)
    exa_free(tr.tree1);
  
  if(tr.constraintVector != NULL)
    exa_free(tr.constraintVector); 
  if(tr.nodeBaseAddress != NULL)
    exa_free(tr.nodeBaseAddress);


  // free parsimony related stuff
  if(parsimonyEnabled)
    {
      exa_free(tr.parsimonyScore);
      for(int i = 0; i < getNumberOfPartitions(); ++i)
	{
	  pInfo *partition = getPartition(i); 
	  exa_free(partition->parsVect);
	}
      exa_free(tr.ti); 
    }


  exa_free(tr.td[0].ti);    
  exa_free(tr.td[0].executeModel);    
  exa_free(tr.td[0].parameterValues);    

  exa_free(tr.nodep); 

  int numPart = getNumberOfPartitions();
  for(int i = 0; i < numPart ;++i)
    {
      pInfo *partition = getPartition(i); 

      exa_free(partition->frequencyGrouping); 
      exa_free(partition->symmetryVector);
      exa_free(partition->frequencies);
      exa_free(partition->empiricalFrequencies); 
      exa_free(partition->tipVector); 
      exa_free(partition->perSiteRates); 
      exa_free(partition->gammaRates); 
      // TODO xvector 
      exa_free(partition->yVector);
      exa_free(partition->xSpaceVector); 
      exa_free(partition->sumBuffer); 

      exa_free(partition->wgt); 
      exa_free(partition->rateCategory); 
      exa_free(partition->partitionName); 
      exa_free(partition->EIGN); 
      exa_free(partition->left); 
      exa_free(partition->right); 
      exa_free(partition->EI); 
      exa_free(partition->EV); 
      exa_free(partition->globalScaler);       
      exa_free(partition->substRates);       
#if HAVE_PLL != 0      
      exa_free(partition->ancestralBuffer); 
      exa_free(partition); 
#endif
    }

#if HAVE_PLL != 0 
  exa_free(partitions.partitionData); 
#else 
  exa_free(tr.executeModel) ;  
  exa_free(tr.partitionContributions); 
  exa_free(tr.fracchanges); 

#endif

}


void TreeAln::initializeFromByteFile(std::string _byteFileName)
{
#if HAVE_PLL != 0
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
  
  initializeTree(&tr, &adef);   
#endif  

  for(int i = 0; i < getNumberOfPartitions(); ++i)
    {
      // TODO empirical freqs? 

      // TODO aa? 

      pInfo *partition =  getPartition(i);
      assert(partition->dataType == DNA_DATA);

      setFrequencies(std::vector<double>(partition->states, 1.0 / (double)partition->states ), i); 
      int num = numStateToNumInTriangleMatrix(partition->states); 
      setRevMat(std::vector<double>( num , 1.0 ), i); 
    }
}


double TreeAln::getTreeLengthExpensive() const
{
  std::vector<Branch> branches = extractBranches(); 
  
  assert(getNumBranches() == 1 ); 

  double result = 1; 
  for(auto b : branches)
    result *= b.getLength(); 

  return result; 
}



void TreeAln::extractHelper( nodeptr p , std::vector<Branch> &result, bool isStart) const 
{
  Branch b = Branch(p->number, p->back->number,  p->back->z[0]);       
  result.push_back(b);

  if(not isStart && isTipNode(p))
    return; 

  assert(getNumBranches() == 1 ); 
  
  for(nodeptr q =  p->next; p != q ; q = q->next)        
    extractHelper( q->back, result, false);    
}


std::vector<Branch> TreeAln::extractBranches() const 
{
  std::vector<Branch> result; 
  extractHelper(tr.nodep[1]->back, result, true);
  return result; 
}


// void TreeAln::verifyTreeLength() const
// {
// #if TODO 
//   double tlVerified = getTreeLengthExpensive() ;
//   assert(treeLength == tlVerified); 
// #endif
// }

void TreeAln::enableParsimony()
{  
#if HAVE_PLL == 0
  // assert the parsimony x-ints are set correctly 
  for(nat i = 1; i < getNumberOfNodes() + 1 ; ++i)
    {
      tr.nodep[i]->xPars = 1; 
      tr.nodep[i]->next->xPars = 0; 
      tr.nodep[i]->next->next->xPars = 0; 
    }
  
  // allocateParsimonyDataStructures(tr);   
#else 
  allocateParsimonyDataStructures(&tr, &partitions);   
#endif
}





/** 
    @brief standard values for the tree 
*/ 
void TreeAln::initDefault()
{   
  tr.likelihood =  0 ; 
  tr.doCutoff = TRUE;
  tr.secondaryStructureModel = SEC_16; /* default setting */
  tr.searchConvergenceCriterion = FALSE;
  tr.rateHetModel = GAMMA; 
  tr.multiStateModel  = GTR_MULTI_STATE;
#if HAVE_PLL == 0 
  tr.useGappedImplementation = FALSE;
  tr.saveBestTrees          = 0;
#endif
  tr.saveMemory = FALSE;
  tr.manyPartitions = FALSE;
  tr.categories             = 25;
  tr.grouped = FALSE;
  tr.constrained = FALSE;
  tr.gapyness               = 0.0; 
  tr.useMedian = FALSE;
  tr.mxtips = 0; 
}


void TreeAln::clipNodeDefault(nodeptr p, nodeptr q )
{
  clipNode(p,q, TreeAln::initBL);
}


void TreeAln::clipNode(nodeptr p, nodeptr q, double z)
{
  p->back = q; 
  q->back = p; 
  assert(getNumBranches() == 1);     
  p->z[0] = p->back->z[0]  = z; 
}


void TreeAln::unlinkTree()
{
#ifdef DEBUG_LINK_INFO
  cout << "unlinking everything" << endl; 
#endif

  // unlink tips 
  for(int i = 1; i < tr.mxtips+1; ++i)
    {
      nodeptr p = tr.nodep[i]; 
      p->back = NULL; 
    }

  for(int i = tr.mxtips + 1; i < 2 * tr.mxtips; ++i)
    {      
      nodeptr p = tr.nodep[i]; 
      p->back = NULL; 
      p->next->back = NULL; 
      p->next->next->back = NULL; 
    }

}


nodeptr TreeAln::getUnhookedNode(int number)
{  
  nodeptr p = tr.nodep[number]; 

  if(isTip(number, tr.mxtips) )
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

      for(int j = 0; j < partitionRhs->states; ++j )
	partitionLhs->frequencies[j] = partitionRhs->frequencies[j]; 
      
      for(nat j = 0 ; j < numStateToNumInTriangleMatrix(partitionRhs->states); ++j)
	partitionLhs->substRates[j] = partitionRhs->substRates[j]; 

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
  for(int i = mxtips+1 ; i < 2* tr.mxtips-1; ++i)
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

  for(int i = mxtips+1; i < 2 * tr.mxtips-1; ++i)
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

  tr.start = tr.nodep[rhsTree->start->number];   
  return *this; 
}



int TreeAln::getNumBranches() const
{
#if HAVE_PLL != 0
  return partitions.perGeneBranchLengths  ? partitions.numberOfPartitions : 1 ; 
#else 
  return tr.numBranches; 
#endif
}



int TreeAln::getNumberOfPartitions() const
{
#if HAVE_PLL != 0
  return partitions.numberOfPartitions; 
#else 
  return tr.NumberOfModels; 
#endif
}



// usefull stuff 
pInfo* TreeAln::getPartition(int model)  const
{
  assert(model < getNumberOfPartitions()); 
  
#if HAVE_PLL != 0  
  return partitions.partitionData[model]; 
#else 
  return tr.partitionData + model; 
#endif
}




// this is BS 
// int TreeAln::accessExecModel(int model) const
// {
// #if HAVE_PLL != 0
//   return partitions.partitionData[model]->executeModel; 
// #else 
//   return tr.executeModel[model]; 
// #endif
// }


// double TreeAln::accessPartitionLH(int model) const 
// {
// #if HAVE_PLL != 0
//   return partitions.partitionData[model]->partitionLH; 
// #else  
//   return tr.perPartitionLH[model]; 
// #endif
// }



std::vector<bool> TreeAln::getExecModel() const 
{ 
  std::vector<bool> result; 
  for(int i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0  
      result.push_back(partitions.partitionData[i]->executeModel); 
#else  
      result.push_back(tr.executeModel[i]); 
#endif
    }
   return result; 
}
 
void TreeAln::setExecModel(const std::vector<bool>  &modelInfo)
{
  for(int i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0 
      partitions.partitionData[i]->executeModel =  modelInfo[i] ?  true : false; 
#else 
      tr.executeModel[i] = modelInfo[i]; 
#endif
    }
}

std::vector<double> TreeAln::getPartitionLnls() const
{
  std::vector<double> result; 
  for(int i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0 
      result.push_back(partitions.partitionData[i]->partitionLH); 
#else 
      result.push_back(tr->perPartitionLH[i]); 
#endif
    }
  return result; 
}
 
void TreeAln::setPartitionLnls(const std::vector<double> partitionLnls)  
{
  for(int i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0
      partitions.partitionData[i]->partitionLH = partitionLnls[i]; 
#else 
      tr.perPartitionLH[i];  = partitionLnls[i]; 
#endif
    }
} 


void TreeAln::initRevMat(int model)
{
#if HAVE_PLL != 0
  initReversibleGTR(getTr(), getPartitionsPtr() , model); 
#else 
  initReversibleGTR(getTr(), model); 
#endif
}


void TreeAln::setFrequencies(const std::vector<double> &values, int model)
{
  // tout << "setting frequencies "; 
  // for_each(values.begin(), values.end(), [](double d) {tout << std::setprecision(3) << d << "," ;  } );
  // tout << std::endl; 

  assert( BoundsChecker::checkFrequencies(values) ) ;    
  auto partition = getPartition(model); 
  memcpy( partition->frequencies, &(values[0]), partition->states * sizeof(double)); 
  initRevMat(model);   
}


void TreeAln::setRevMat(const std::vector<double> &values, int model)
{
  // std::cout << "checking "; for_each(values.begin(), values.end(), [](double d){std::cout << d << "," ;}) ; std::cout  << std::endl; 
  assert( BoundsChecker::checkRevmat(values) ); 
  auto partition = getPartition(model) ; 
  memcpy(partition->substRates, &(values[0]), numStateToNumInTriangleMatrix(  partition->states) * sizeof(double)); 
  initRevMat(model); 
}



void TreeAln::setBranch(const Branch& branch)
{
  assert(BoundsChecker::checkBranch(branch)); 
  auto p = branch.findNodePtr(*this);
  p->z[0] = p->back->z[0] = branch.getLength();
}


void TreeAln::setAlpha(double alpha,  int model)
{
  assert(BoundsChecker::checkAlpha(alpha)); 
  auto p = getPartition(model); 
  p->alpha = alpha; 
  discretizeGamma(model); 
}


/**
   @brief makes the discrete categories for the gamma
   distribution. Has to be called, if alpha was updated.   

 */ 
void TreeAln::discretizeGamma(int model)
{
  pInfo *partition =  getPartition(model); 
  makeGammaCats(partition->alpha, partition->gammaRates, 4, tr.useMedian);
}




// initialization code from the PLL in future version, we could also
// support plain files again. Just do not really know what to do with
// it =/
#if HAVE_PLL  != 0 

void TreeAln::initializeTreePLL(std::string byteFileName)
{
  partitionList *pl = getPartitionsPtr();
  pl->partitionData = (pInfo**)exa_malloc(NUM_BRANCHES*sizeof(pInfo*));

  double **empFreq = NULL; 
  bool perGeneBL = false; 
  initializePartitionsPLL(byteFileName, &empFreq, perGeneBL );

  initializePartitionsSequential(getTr(), getPartitionsPtr());
  initModel(getTr(), empFreq, getPartitionsPtr());
  
  for(int i = 0; i < pl->numberOfPartitions; ++i)
    exa_free(empFreq[i]);
  exa_free(empFreq) ;
} 


void TreeAln::initializePartitionsPLL(std::string byteFileName, double ***empiricalFrequencies, bool multiBranch)
{
  unsigned char *y;

  FILE 
    *byteFile = myfopen(byteFileName.c_str(), "rb");	 

  myBinFread(&(tr.mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr.originalCrunchedLength), sizeof(int), 1, byteFile);
  myBinFread(&(partitions.numberOfPartitions),  sizeof(int), 1, byteFile);
  myBinFread(&(tr.gapyness),            sizeof(double), 1, byteFile);


  partitions.perGeneBranchLengths = multiBranch ? 1 : 0 ;

  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed RAxML-specific */

  tr.aliaswgt                   = (int *)exa_malloc((size_t)tr.originalCrunchedLength * sizeof(int));
  myBinFread(tr.aliaswgt, sizeof(int), tr.originalCrunchedLength, byteFile);	       

  tr.rateCategory    = (int *)    exa_malloc((size_t)tr.originalCrunchedLength * sizeof(int));	  

  tr.patrat          = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double));
  tr.patratStored    = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double)); 
  tr.lhs             = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double)); 

  *empiricalFrequencies = (double **)exa_malloc(sizeof(double *) * partitions.numberOfPartitions);

  y = (unsigned char *)exa_malloc(sizeof(unsigned char) * ((size_t)tr.originalCrunchedLength) * ((size_t)tr.mxtips));
  tr.yVector = (unsigned char **)exa_malloc(sizeof(unsigned char *) * ((size_t)(tr.mxtips + 1)));

  for( nat i = 1; i <= (size_t)tr.mxtips; i++)
    tr.yVector[i] = &y[(i - 1) *  (size_t)tr.originalCrunchedLength]; 

  setupTree(&tr, FALSE, &partitions);

  for(int i = 0; i < partitions.numberOfPartitions; i++)
    partitions.partitionData[i]->executeModel = TRUE;

  /* data structures for convergence criterion need to be initialized after! setupTree */

  for( nat i = 1; i <= (size_t)tr.mxtips; i++)
    {
      int len;
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr.nameList[i] = (char*)exa_malloc(sizeof(char) * (size_t)len);
      myBinFread(tr.nameList[i], sizeof(char), len, byteFile);
      /*printf("%s \n", tr.nameList[i]);*/
    }  

  for( nat i = 1; i <= (size_t)tr.mxtips; i++)
    addword(tr.nameList[i], tr.nameHash, i);

  for(int model = 0; model < partitions.numberOfPartitions; model++)
    {
      int 
	len;

      pInfo 
	*p = partitions.partitionData[model];

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

      (*empiricalFrequencies)[model] = (double *)exa_malloc(sizeof(double) * (size_t)partitions.partitionData[model]->states);
      myBinFread((*empiricalFrequencies)[model], sizeof(double), partitions.partitionData[model]->states, byteFile);
    }

  myBinFread(y, sizeof(unsigned char), ((size_t)tr.originalCrunchedLength) * ((size_t)tr.mxtips), byteFile);

  fclose(byteFile);
}


#endif


std::ostream& operator<< (std::ostream& out,  TreeAln&  traln)
{
  TreePrinter tp(true, false, false); 
  return out << tp.printTree(traln); 
}



void TreeAln::collapseBranch(Branch b)
{
  assert(getNumBranches() == 1 ); 
  nodeptr p = b.findNodePtr( *this); 
  p->z[0] = p->back->z[0] = TreeAln::zZero;   
}


bool TreeAln::isCollapsed(Branch b ) 
{
  assert(getNumBranches() == 1 ); 
  nodeptr p = b.findNodePtr(*this); 
  return p->z[0] >=  BoundsChecker::zMax ; 
}


void TreeAln::setBranchLengthUnsafe(Branch b ) 
{
  assert(getNumBranches() == 1 ); 
  nodeptr p = b.findNodePtr(  *this);   
  p->z[0] = p->back->z[0] = b.getLength();   
} 


std::vector<double> TreeAln::getRevMat(int model) const 
{
  std::vector<double> result; 
  pInfo *partition = getPartition(model); 
  double sum = 0; 
  for(nat i =0 ; i < numStateToNumInTriangleMatrix(partition->states); ++i)
    {
      result.push_back(partition->substRates[i]);
      sum += result[i]; 
    }

  for_each(result.begin(), result.end(), [&](double &v) { v /= sum ; }) ; 

  return result; 
}


std::vector<double> TreeAln::getFrequencies(int model) const
{
  std::vector<double> result; 
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


nat numStateToNumInTriangleMatrix(int numStates)  
{  
  return (  numStates * numStates - numStates) / 2 ; 
}




nodeptr TreeAln::getNode(nat elem) const  
{
  if( not (elem != 0 && elem < getNumberOfNodes() + 2 ))
    { 
      std::cout << "bug: attempted to get node " << elem << std::endl; assert(elem != 0 && elem <= getNumberOfNodes() + 1 ) ;
    } 

  return  getTr()->nodep[elem] ; 
}
