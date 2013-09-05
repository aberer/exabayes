#include <sstream>
#include <cassert>
#include <cstring>
#include <unordered_set>

#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/AbstractParameter.hpp"
#include "TreeRandomizer.hpp"
#include "axml.h"
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"
#include "BoundsChecker.hpp"
#include "TreePrinter.hpp"

const double TreeAln::zZero = BoundsChecker::zMax + ( 1 - BoundsChecker::zMax) / 2 ; 
const double TreeAln::initBL = 0.75;
const double TreeAln::problematicBL = std::numeric_limits<double>::max();

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
      for(nat i = 0; i < getNumberOfPartitions(); ++i)
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

  adef.perGeneBranchLengths   = TRUE;   

  adef.useCheckpoint          = FALSE;
  
  initializeTree(&tr, &adef);   
#endif  


  // set default values (will be overwritten later, if necessary )
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
      pInfo *partition =  getPartition(i);
      assert(partition->dataType == DNA_DATA);

      setFrequencies(std::vector<double>(partition->states, 1.0 / (double)partition->states ), i); 
      int num = numStateToNumInTriangleMatrix(partition->states); 
      setRevMat(std::vector<double>( num , 1.0 ), i); 
      setAlpha(1, i); 
    }

  // HACK: initialize the tree somehow 
  randCtr_t s; 
  s.v[0] = 0 ; 
  s.v[1] = 1 ; 
  Randomness r(s); 
  TreeRandomizer::randomizeTree(*this,r ); 

  // evaluate once, to activate the likelihood arrays  
  auto start = tr.start; 
#if HAVE_PLL != 0  
  evaluateGeneric(&tr, getPartitionsPtr(), start, true); 
#else 
  evaluateGeneric(&tr, start, true); 
#endif   
}


void TreeAln::extractHelper( nodeptr p , std::vector<BranchLength> &result, bool isStart, const AbstractParameter* param ) const 
{
  // TODO inefficient 
#ifdef EFFICIENT
  assert(0); 
#endif
  auto tmp = BranchPlain(p->number, p->back->number); 
  auto b = getBranch(tmp, param);

  result.push_back(b);

  if(not isStart && isTipNode(p))
    return; 

  for(nodeptr q =  p->next; p != q ; q = q->next)        
    extractHelper( q->back, result, false, param);
}

void TreeAln::extractHelper( nodeptr p , std::vector<BranchLengths> &result, bool isStart, const std::vector<AbstractParameter*> &params ) const 
{
  // TODO inefficient 
#ifdef EFFICIENT
  assert(0); 
#endif
  auto tmp = BranchPlain(p->number, p->back->number); 
  auto b = getBranch(tmp, params);

  result.push_back(b);

  if(not isStart && isTipNode(p))
    return; 

  for(nodeptr q =  p->next; p != q ; q = q->next)        
    extractHelper( q->back, result, false, params);
}

void TreeAln::extractHelper( nodeptr p , std::vector<BranchPlain> &result, bool isStart ) const 
{
  // TODO inefficient 
#ifdef EFFICIENT
  assert(0); 
#endif
  auto tmp = BranchPlain(p->number, p->back->number); 

  // auto b = getBranch(tmp);

  result.push_back(tmp);
  if(not isStart && isTipNode(p))
    return; 
  for(nodeptr q =  p->next; p != q ; q = q->next)        
    extractHelper( q->back, result, false);
}


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
  tr.numBranches = getNumberOfPartitions();
#else 
  partitions.perGeneBranchLengths = TRUE; 
#endif

  // TODO 
  tr.manyPartitions = FALSE;
  tr.saveMemory = FALSE;

  tr.categories             = 25;
  tr.grouped = FALSE;
  tr.constrained = FALSE;
  tr.gapyness               = 0.0; 
  tr.useMedian = FALSE;
  tr.mxtips = 0; 
}


void TreeAln::clipNodeDefault(nodeptr p, nodeptr q )
{
  clipNode(p,q);  
  // use a problematic branch length to ensure, that it was overridden
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    p->z[i] = p->back->z[i] = TreeAln::initBL; 
}


void TreeAln::clipNode(nodeptr p, nodeptr q)
{
  p->back = q; 
  q->back = p; 

#ifdef EFFICIENT
  // the following only done for error checking 
  assert(0);
#endif
  
  // would be cooler, if we could keep that 
  // for(nat i = 0; i < getNumberOfPartitions(); ++i )
  //   p->z[i] = p->back->z[i]  = std::numeric_limits<double>::max(); 
}


void TreeAln::detachNode(nodeptr p)
{
  p->next->back = NULL; 
  p->next->next->back = NULL; 
}

void TreeAln::unlinkNode(nodeptr p )
{
  p->back = NULL; 
  detachNode(p); 
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
    unlinkNode(tr.nodep[i]); 
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


bool TreeAln::operator==(const TreeAln& rhs)
{
#if 0 
  // notice: this is untested! 
  assert(0); 

  bool result = true; 

  // create some artificial parameters 
  std::vector<BranchLengthsParameter> realParams; 
  std::vector<AbstractParameter*> params; 
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
      realParams.emplace_back(i, i);
      auto &param  =  *(realParams.rbegin()); 
      param.addPartition(i);
    }
  for(auto & p: realParams)
    params.push_back(&p); 


  auto rBranches = rhs.extractBranches(params); 
  auto lBranches = extractBranches(params); 

  auto equals = [](const Branch& a, const Branch &b ) { return a.equals(b, BranchEqualFlag::WITH_LENGTH) ;}; 
  auto theHash = [](const Branch &a) { return std::hash<nat>()(a.getPrimNode()) ^ std::hash<nat>()(a.getSecNode()); }; 
  std::unordered_set<Branch, decltype(theHash), decltype(equals)> lSet(lBranches.size(), theHash, equals); 
  std::unordered_set<Branch, decltype(theHash), decltype(equals)> rSet(rBranches.size(), theHash, equals); 
  for(auto &b : rBranches)
    rSet.insert(b); 
  for(auto &b : lBranches)
    lSet.insert(b); 
  
  for(auto &b : lBranches)
    result &= rSet.find(b)!= rSet.end(); 
  result &= lBranches.size() == rBranches.size(); 

  
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
      const auto rp =  rhs.getPartition(i); 
      const auto lp = getPartition(i); 

      result &= rp->dataType == lp->dataType; 

      if(result)
	{
	  // check frequencies 
	  for(int j = 0; j < lp->states; ++j)
	    result &= rp->frequencies[j] == lp->frequencies[j]; 
      
	  // check revmat 
	  for(nat j = 0; j < numStateToNumInTriangleMatrix(lp->states) ; ++j) 
	    result &= lp->substRates[j] == rp->substRates[j]; 
	}

      // check alpha
      result &= rp->alpha == lp->alpha; 

    }

  return result; 
#else 
  assert(0); 
  return false; 
#endif
}



#define UNSAFE_EXACT_TREE_COPY

void TreeAln::copyModel(const TreeAln& rhs)  
{  
  // assert(&rhs != this); 
  if(&rhs == this)
    return; 

  // copy partition parameters 
  for(nat i = 0; i < rhs.getNumberOfPartitions(); ++i)
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
  auto *rhsTree = rhs.getTr(); 
  auto *thisTree = getTr();


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
}



nat TreeAln::getNumBranches() const
{
#if HAVE_PLL != 0
  return partitions.perGeneBranchLengths  ? partitions.numberOfPartitions : 1 ; 
#else 
  return tr.numBranches; 
#endif
}


nat TreeAln::getNumberOfPartitions() const
{
#if HAVE_PLL != 0
  return partitions.numberOfPartitions; 
#else 
  return tr.NumberOfModels; 
#endif
}


pInfo* TreeAln::getPartition(nat model)  const
{
  assert(model < getNumberOfPartitions()); 
  
#if HAVE_PLL != 0  
  return partitions.partitionData[model]; 
#else 
  return tr.partitionData + model; 
#endif
}


std::vector<bool> TreeAln::getExecModel() const 
{ 
  std::vector<bool> result; 
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
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
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
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
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0 
      result.push_back(partitions.partitionData[i]->partitionLH); 
#else 
      result.push_back(tr.perPartitionLH[i]); 
#endif
    }
  return result; 
}
 
void TreeAln::setPartitionLnls(const std::vector<double> partitionLnls)  
{
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0
      partitions.partitionData[i]->partitionLH = partitionLnls[i]; 
#else 
      tr.perPartitionLH[i]  = partitionLnls[i]; 
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

  bool valuesOkay = BoundsChecker::checkRevmat(values); 
  if(not valuesOkay)
    {
      tout << "Problem with substitution parameter: "  << MAX_SCI_PRECISION << values << std::endl;  
      assert( valuesOkay ); 
    }

  auto partition = getPartition(model) ; 
  memcpy(partition->substRates, &(values[0]), numStateToNumInTriangleMatrix(  partition->states) * sizeof(double)); 
  initRevMat(model); 
}


void TreeAln::setBranch(const BranchLength& branch, const AbstractParameter* param)
{
  assert(BoundsChecker::checkBranch(branch)); 
  assert(branch.exists(*this)); 

  auto p = branch.findNodePtr(*this);
  for(auto &partition : param->getPartitions() )
    {
      double length = branch.getLength(); 
      p->z[partition] = p->back->z[partition] = length; 
    }
}


void TreeAln::setBranch(const BranchLengths& branch, const std::vector<AbstractParameter*>params)
{
  assert(BoundsChecker::checkBranch(branch)); 
  assert(branch.exists(*this)); 

  auto p = branch.findNodePtr(*this); 
  
  for(auto &param : params)
    {
      double length = branch.getLengths().at(param->getIdOfMyKind());
      for(auto &partition : param->getPartitions())
	p->z[partition] = p->back->z[partition] = length; 
    }
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
  // partitions.perGeneBranchLengths = TRUE; 
  partitionList *pl = getPartitionsPtr();
  pl->partitionData = (pInfo**)exa_malloc(NUM_BRANCHES*sizeof(pInfo*));

  double **empFreq = NULL; 
  bool perGeneBL = true; 
  initializePartitionsPLL(byteFileName, &empFreq, perGeneBL );

  initializePartitionsSequential(getTr(), getPartitionsPtr());
  initModel(getTr(), empFreq, getPartitionsPtr());
  
  for(int i = 0; i < pl->numberOfPartitions; ++i)
    exa_free(empFreq[i]);
  exa_free(empFreq) ;
  partitions.perGeneBranchLengths = TRUE; 
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


std::ostream& operator<< (std::ostream& out,  const TreeAln&  traln)
{
  TreePrinter tp(true, false, false); 
  return out << tp.printTree(traln); 
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


std::pair<BranchPlain,BranchPlain> TreeAln::getDescendents(const BranchPlain &b) const
{
  auto p = b.findNodePtr(*this); 
  return 
    std::make_pair ( BranchPlain(p->next->number, p->next->back->number), 
		     BranchPlain(p->next->next->number, p->next->next->back->number) ); 
} 


 
// for debug 
void TreeAln::printArrayStart(std::ostream &out, nat length )
{
  out << SOME_SCI_PRECISION; 
  auto numTax = getNumberOfTaxa(); 
  for(nat i =   0  ; i < getNumberOfInnerNodes( )   ; ++i)
    {
      out << "[ array " << i + numTax  <<  " ]"; 
	
      for(nat ctr = 0 ; ctr < length  ; ++ctr)
	{
	  for(nat j = 0 ; j < 4; ++j)
	    {
	      auto partition = getPartition(0); 
	      out << partition->xVector[i][ctr * 4 + j] << "," ; 
	    }
	  out << "\t" ;
	}
      out << std::endl; 
    } 
}




std::ostream& operator<<(std::ostream& out, pInfo& rhs)
{
  assert(rhs.dataType == DNA_DATA); 
  out << "{DNA,pi=(" ; 
  std::copy(rhs.frequencies, rhs.frequencies + 4, std::ostream_iterator<double>(tout,",")); 
  out << "),r=("; 
  std::copy(rhs.substRates, rhs.substRates+ 6, std::ostream_iterator<double> (tout, ",")); 
  out << "),alpha=" << rhs.alpha << "}"; 
  return out ; 
}


double TreeAln::getMeanSubstitutionRate(const std::vector<nat> &partitions) const 
{
  double result = 0; 
  // TODO  
  // tout << "not sure about this" << std::endl; 

#if HAVE_PLL == 0
  for(auto &p :partitions)
    result += tr.fracchanges[p] * tr.partitionContributions[p];   
#else 
  for(auto &p : partitions)
    {
      auto partition =  getPartition(p);
      result += partition->fracchange * partition->partitionContribution;  
    }
#endif
  // tout << "fracchange " << result << " for partitions " << partitions << std::endl; 
  return result; 
}



// Branch TreeAln::getBranch(const Branch& b, AbstractParameter* const &param) const
// {
//   return getBranch(b.findNodePtr(*this) , param);
// }


// Branch TreeAln::getBranch(const Branch &b , const std::vector<AbstractParameter*> &params) const 
// {
//   return getBranch(b.findNodePtr(*this),params);
// }


// // Branch TreeAln::getBranch(nodeptr p, AbstractParameter* const &param) const
// // {
// //   nat id = param->getIdOfMyKind();
// //   // std::vector<double> l(id+1, TreeAln::problematicBL); // =/ 
  
// //   auto l = 

// //   nat aPartition = param->getPartitions()[0];
// //   l.at(id) = p->z[aPartition];
// //   return Branch(p->number,p->back->number, l);
// // }


// Branch TreeAln::getBranch(nodeptr p, AbstractParameter* const &param) const 
// {
//   auto params = std::vector<AbstractParameter*> {param}; 
//   return getBranch(p, params); 
// }


// Branch TreeAln::getBranch(nodeptr p, const std::vector<AbstractParameter*> &params) const 
// {
//   nat highest = 0; 
//   for(auto &p : params )
//     {
//       nat elem = p->getIdOfMyKind(); 
//       if(highest < elem )
// 	highest = elem; 
//     }

//   auto l = std::vector<double>{}; 
//   for(auto &param: params)
//     {
//       nat aPartition = param->getPartitions()[0];
//       nat id = param->getIdOfMyKind(); 
//       if(id >= l.size())
// 	l.resize(id + 1); 
//       l.at(id) = p->z[aPartition];
//     } 
//   return Branch(p->number, p->back->number, l); 
// }



std::vector<std::string> TreeAln::getNameMap() const
{
  auto result = std::vector<std::string>{};
  auto tr = getTr();   

  result.push_back("error"); 
  for(nat i = 1 ; i < getNumberOfTaxa()+1; ++i)
    result.push_back(tr->nameList[i]);

  return result; 
} 


nat TreeAln::getDepth(const BranchPlain &b) const 
{
  if(b.isTipBranch(*this))
    return 1; 
  else 
    {
      auto desc = getDescendents(b.getInverted());
      return std::max(getDepth(desc.first), getDepth(desc.second)) + 1; 
    }
}


std::vector<nat> TreeAln::getLongestPathBelowBranch(const BranchPlain &b) const 
{
  if(b.isTipBranch(*this))
    return std::vector<nat>{b.getSecNode()}; 
  else 
    {
      auto desc = getDescendents(b.getInverted());
      auto resA = getLongestPathBelowBranch(desc.first);
      auto resB = getLongestPathBelowBranch(desc.second);
      if(resA.size() > resB.size())
	{
	  resA.push_back(b.getSecNode()); 
	  return resA; 
	}
      else 
	{
	  resB.push_back(b.getSecNode()); 
	  return resB; 
	}
    }
}


std::vector<nat> TreeAln::getNeighborsOfNode( nat node ) const 
{
  auto result = std::vector<nat>{};
  auto nodePtr = getNode(node); 
  if(isTipNode(nodePtr))
    result.push_back(nodePtr->back->number);
  else 
    {
      result = 
	{
	  nat(nodePtr->back->number), 
	  nat(nodePtr->next->back->number), 
	  nat(nodePtr->next->next->back->number)
	} ; 
    }
  
  return result; 
} 


std::vector<BranchPlain> TreeAln::extractBranches() const 
{
  auto result = std::vector<BranchPlain>{}; 
  extractHelper(tr.nodep[1]->back, result, true);
  return result;
}


std::vector<BranchLength> TreeAln::extractBranches(const AbstractParameter* param) const 
{
  auto result = std::vector<BranchLength>{}; 
  extractHelper(tr.nodep[1]->back, result, true, param);
  return result;
}


std::vector<BranchLengths> TreeAln::extractBranches(const std::vector<AbstractParameter*> &params) const 
{
  auto result = std::vector<BranchLengths>{}; 
  extractHelper(tr.nodep[1]->back, result, true, params);
  return result;
}


BranchPlain TreeAln::getAnyBranch() const 
{
  return BranchPlain(tr.nodep[1]->number, tr.nodep[1]->back->number); 
} 

BranchLength TreeAln::getBranch(const nodeptr &p,  const AbstractParameter *param) const
{
  return getBranch(BranchPlain(p->number,p->back->number), param); 
}


BranchLengths TreeAln::getBranch(const nodeptr& p, const std::vector<AbstractParameter*> &params) const 
{
  return getBranch(BranchPlain(p->number, p->back->number), params); 
}


BranchLength TreeAln::getBranch(const BranchPlain &branch,  const AbstractParameter *param) const
{
  auto result = branch.toBlDummy();
  result.extractLength(*this, branch, param); 
  return result; 
}


BranchLengths TreeAln::getBranch(const BranchPlain& branch, const std::vector<AbstractParameter*> &params) const 
{
  auto result = branch.toBlsDummy();
  result.extractLength(*this, branch, params);
  return result; 
}

