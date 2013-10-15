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
#include "TreeInitializer.hpp"

const double TreeAln::initBL = 0.1;


TreeAln::TreeAln()
  : parsimonyEnabled(true)
  , mode(RunModes::NOTHING)
{
  memset(&tr,0,sizeof(tree)); 

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

  tr.manyPartitions = FALSE;
  tr.saveMemory = FALSE;
  tr.categories             = 25;
  tr.grouped = FALSE;
  tr.constrained = FALSE;
  tr.gapyness               = 0.0; 
  tr.useMedian = FALSE;
  tr.mxtips = 0; 
}


void TreeAln::initializeNumBranchRelated()
{
  NUM_BRANCHES = getNumberOfPartitions();
  tr.zqr = new double[NUM_BRANCHES]; 
  tr.currentZQR = new double[NUM_BRANCHES]; 
  tr.currentLZR = new double[NUM_BRANCHES];
  tr.currentLZQ = new double[NUM_BRANCHES];
  tr.currentLZS = new double[NUM_BRANCHES];
  tr.currentLZI = new double[NUM_BRANCHES];
  tr.lzs = new double[NUM_BRANCHES];
  tr.lzq = new double[NUM_BRANCHES];
  tr.lzr = new double[NUM_BRANCHES];
  tr.lzi = new double[NUM_BRANCHES];
  tr.coreLZ = new double[NUM_BRANCHES]; 
  tr.curvatOK = new boolean[NUM_BRANCHES];  

  tr.partitionSmoothed = new boolean[NUM_BRANCHES];
  tr.partitionConverged = new boolean[NUM_BRANCHES]; 

  nat numTax = getNumberOfTaxa() ; 

  for(nat i = 0; i < numTax ; ++i)
    {
      tr.td[0].ti[i].qz  = new double[NUM_BRANCHES]; 
      tr.td[0].ti[i].rz  = new double[NUM_BRANCHES]; 
    }

  for(nat i = 0; i < numTax + 3 * (numTax - 1)  ; ++i )
    {
      auto node = tr.nodeBaseAddress + i; 
      node->z = new double[NUM_BRANCHES]; 
    }
}


void TreeAln::clearMemory()
{
  for(nat i = 0; i < getNumberOfPartitions() ; ++i)
    {
      auto partition = getPartition(i);
      for(nat j = 0; j < getNumberOfInnerNodes() ; ++j)
	{
	  if(partition->xSpaceVector[j] != 0 )
	    {
	      // tout << "freeing " << j <<   " which had length "<<partition->xSpaceVector[j] << std::endl; 
	      exa_free( partition->xVector[j] );
	      partition->xVector[j] = NULL; 
	      partition->xSpaceVector[j] = 0; 
	    }
	}
    }
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

  clearMemory(); 

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




void TreeAln::initializeFromByteFile(std::string _byteFileName, RunModes flags)
{
  mode = flags; 

  if( ( flags & RunModes::PARTITION_DISTRIBUTION)  !=  RunModes::NOTHING)
    tr.manyPartitions = TRUE; 
  
  if( (flags & RunModes::MEMORY_SEV) != RunModes::NOTHING)
    tr.saveMemory = TRUE; 
  
  auto ti = TreeInitializer{};  
  ti.unifiedInitializePartitions(*this, _byteFileName ); 
  
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

  initializeNumBranchRelated();
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

#else 
  allocateParsimonyDataStructures(&tr, &partitions);   
#endif
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
  
  assert(not isTipNode(p)); 

  auto pn = p->next,
    pnn = p->next->next; 

  return 
    std::make_pair ( BranchPlain(pn->number, pn->back->number), 
		     BranchPlain(pnn->number, pnn->back->number) ); 
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

  return result; 
}


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


BranchPlain TreeAln::getAnyBranch() const 
{
  // assert(0); 
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

nat TreeAln::getNumberOfAssignedSites(nat model) const 
{
  auto partition =  getPartition (model); 
#if HAVE_PLL != 0
  nat length = partition->upper - partition->lower; 
#else 
  nat length = partition->width; 
#endif

  return length; 
}


std::vector<BranchPlain> TreeAln::extractBranches() const
{
  auto result = std::vector<BranchPlain>(); 
  result.reserve(getNumberOfBranches()); 

  int last = int(getNumberOfNodes() + 1); 
  for(int i = getNumberOfTaxa() + 1 ; i < last ; ++i)
    {
      auto node = getNode(i); 
      
      if( node->back->number < i )
	result.emplace_back(i,node->back->number); 
      if(node->next->back->number < i)
	result.emplace_back(i,node->next->back->number); 
      if(node->next->next->back->number < i)
	result.emplace_back(i,node->next->next->back->number); 
    }

  assert(result.size() == getNumberOfBranches());

  return result;
}



std::vector<BranchLength> TreeAln::extractBranches( const AbstractParameter* param ) const 
{
  auto result = std::vector<BranchLength>(); 
  result.reserve(getNumberOfBranches());

  for(int i = getNumberOfTaxa() + 1 ; i < int(getNumberOfNodes() + 1 )   ; ++i)
    {
      auto node = getNode(i); 
      
      if( node->back->number < i )
	result.emplace_back( getBranch(BranchPlain(   i,node->back->number ),param) ); 
      if(node->next->back->number < i)
	result.emplace_back( getBranch(BranchPlain(i,node->next->back->number), param) ); 
      if(node->next->next->back->number < i)
	result.emplace_back( getBranch(BranchPlain( i,node->next->next->back->number ),param) );  
   }

  assert(result.size() == getNumberOfBranches());
  return result;
}

std::vector<BranchLengths> TreeAln::extractBranches( const std::vector<AbstractParameter*> &param ) const 
{
  auto result = std::vector<BranchLengths>(); 
  result.reserve(getNumberOfBranches()); 

  for(int i = getNumberOfTaxa() + 1 ; i < int(getNumberOfNodes() + 1 )  ; ++i)
    {
      auto node = getNode(i); 
      
      if( node->back->number < i )
	result.emplace_back(  getBranch(BranchPlain(   i,node->back->number ),param) ); 
      if(node->next->back->number < i)
	result.emplace_back(  getBranch(BranchPlain(i,node->next->back->number), param)); 
      if(node->next->next->back->number < i)
	result.emplace_back( getBranch(BranchPlain( i,node->next->next->back->number ),param) ); 
    }

  assert(result.size() == getNumberOfBranches());
  return result;
}


std::vector<BranchPlain> TreeAln::getBranchesByDistance(const BranchPlain& branch, nat distance, bool bothSides ) const 
{
  if(distance == 0)
    return { branch };

  auto toCheck = std::vector<BranchPlain>{}; 
  
  if(not isTipNode(branch.findNodePtr(*this))) 
    {
      auto desc = getDescendents(branch); 
      toCheck.push_back(desc.first.getInverted()); 
      toCheck.push_back(desc.second.getInverted()); 
    }

  if(bothSides
     && not isTipNode(branch.getInverted().findNodePtr(*this)))
    {
      auto desc2 = getDescendents(branch.getInverted()); 
      toCheck.push_back(desc2.first.getInverted()); 
      toCheck.push_back(desc2.second.getInverted()); 
    }

  auto myResult = std::vector<BranchPlain>{}; 
  for(auto b : toCheck)
    {
      auto result = getBranchesByDistance(b, distance-1, false); 
      myResult.insert(end(myResult), begin(result), end(result)); 
    }

  return myResult; 
}
