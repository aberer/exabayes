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


void TreeAln::createStandardTree(nat numTax)
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

  tr.manyPartitions = FALSE;
  tr.saveMemory = FALSE;
  tr.categories             = 25;
  tr.grouped = FALSE;
  tr.constrained = FALSE;
  tr.gapyness               = 0.0; 
  tr.useMedian = FALSE;
  tr.mxtips = numTax; 

  setNumberOfPartitions(0);

  // auto tinit = TreeInitializer{}; 
  TreeInitializer::setupTheTree(tr);
}

TreeAln::TreeAln(nat numTax)
  : mode(RunModes::NOTHING)
{
  memset(&tr,0,sizeof(tree)); 
  createStandardTree(numTax);
}


TreeAln::TreeAln(const TreeAln& rhs)
  : mode(rhs.mode)
  , _taxa(rhs._taxa)
{
  tout << "copying a tree "<< std::endl; 
  memset(&tr,0,sizeof(tree));
  createStandardTree(rhs.getNumberOfTaxa());
  
  // copy the topology 
  for(nat i = 1; i < getNumberOfNodes() ; ++i)
    {
      auto rhsP = rhs.getNode(i);
      auto rhsQ = rhsP->back; 
      
      auto p = getUnhookedNode(i); 
      auto q = getUnhookedNode(rhsQ->number);
      clipNode(p,q);

      if(i > getNumberOfTaxa() )
	{
	  p = getUnhookedNode(i);
	  rhsQ = rhsP->next->back; 
	  p = getUnhookedNode(i); 
	  q = getUnhookedNode(rhsQ->number);
	  clipNode(p,q);

	  p = getUnhookedNode(i);
	  rhsQ = rhsP->next->next->back; 
	  p = getUnhookedNode(i); 
	  q = getUnhookedNode(rhsQ->number);
	  clipNode(p,q);
	}
    }

}

void swap(TreeAln &lhs, TreeAln& rhs)
{
  using std::swap; 
  // TODO
  assert(0);
}

TreeAln& TreeAln::operator=(TreeAln rhs)
{
  tout << "using operator= of traln" << std::endl; 
  swap(*this, rhs);
  return *this; 
}

void TreeAln::clearMemory()
{
  for(nat i = 0; i < getNumberOfPartitions() ; ++i)
    {
      auto& partition = getPartition(i);
      for(nat j = 0; j < getNumberOfInnerNodes() ; ++j)
	{
	  if(partition.xSpaceVector[j] != 0 )
	    {
	      // tout << "freeing " << j <<   " which had length "<<partition.xSpaceVector[j] << std::endl; 
	      exa_free( partition.xVector[j] );
	      partition.xVector[j] = NULL; 
	      partition.xSpaceVector[j] = 0; 
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
  exa_free(tr.parsimonyScore);
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
      auto& partition = getPartition(i); 
      if(partition.parsVect != NULL)
	exa_free(partition.parsVect);
    }
  exa_free(tr.ti); 

  exa_free(tr.td[0].ti);    
  exa_free(tr.td[0].executeModel);    
  exa_free(tr.td[0].parameterValues);    

  exa_free(tr.nodep); 

  clearMemory(); 

  int numPart = getNumberOfPartitions();
  for(int i = 0; i < numPart ;++i)
    {
      auto& partition = getPartition(i); 

      exa_free(partition.frequencyGrouping); 
      exa_free(partition.symmetryVector);
      exa_free(partition.frequencies);
      exa_free(partition.empiricalFrequencies); 
      exa_free(partition.tipVector); 
      exa_free(partition.perSiteRates); 
      exa_free(partition.gammaRates); 

      exa_free(partition.yVector);
      exa_free(partition.xSpaceVector); 
      exa_free(partition.sumBuffer); 

      exa_free(partition.wgt); 
      exa_free(partition.rateCategory); 
      exa_free(partition.partitionName); 
      exa_free(partition.EIGN); 
      exa_free(partition.left); 
      exa_free(partition.right); 
      exa_free(partition.EI); 
      exa_free(partition.EV); 
      exa_free(partition.globalScaler);       
      exa_free(partition.substRates);       
    }

#if HAVE_PLL != 0 
  if(getNumberOfPartitions() > 0)
    exa_free(partitions.partitionData); 
#else 
  exa_free(tr.executeModel) ;  
  exa_free(tr.partitionContributions); 
  exa_free(tr.fracchanges); 
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

  
  std::cerr << "Error: did not find unlinked node for " << number << std::endl; 
  assert(0);
  return NULL;
}




#define UNSAFE_EXACT_TREE_COPY

void TreeAln::copyModel(const TreeAln& rhs)  
{  
  // assert(0);
  assert(&rhs != this); 
  // if(&rhs == this)
  //   return; 

  // copy partition parameters 
  for(nat i = 0; i < rhs.getNumberOfPartitions(); ++i)
    {
      auto& partitionRhs = rhs.getPartition(i); 
      auto& partitionLhs = this->getPartition(i); 

      for(int j = 0; j < partitionRhs.states; ++j )
	partitionLhs.frequencies[j] = partitionRhs.frequencies[j]; 
      
      for(nat j = 0 ; j < numStateToNumInTriangleMatrix(partitionRhs.states); ++j)
	partitionLhs.substRates[j] = partitionRhs.substRates[j]; 

      partitionLhs.protModels = partitionRhs.protModels; 	

      partitionLhs.alpha = partitionRhs.alpha; 
      initRevMat(i);
      discretizeGamma(i);       
    }

  this->unlinkTree();
  int mxtips = rhs.getTrHandle().mxtips ;
  auto &rhsTree = rhs.getTrHandle(); 
  auto &thisTree = getTrHandle();

#ifdef UNSAFE_EXACT_TREE_COPY
  // if this works, it is one of the most hackney things, I've ever done... 
  for(int i = mxtips+1 ; i < 2* tr.mxtips-1; ++i)
    {
      nodeptr rhsNode = rhsTree.nodep[i],
	lhsNode = thisTree.nodep[i];       
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
      nodeptr a = rhsTree.nodep[i],
	b = rhsTree.nodep[i]->back; 

#ifdef DEBUG_LINK_INFO
      cout << "hooking up tip " << a->number << " and " << b->number << endl; 
#endif
      
      hookup(this->getUnhookedNode(a->number), this->getUnhookedNode(b->number),a->z, getNumBranches()); 
    }

  for(int i = mxtips+1; i < 2 * tr.mxtips-1; ++i)
    {      
      nodeptr pRhs = rhsTree.nodep[i], 
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


pInfo& TreeAln::getPartition(nat model)  const
{
  assert(model < getNumberOfPartitions()); 
  
#if HAVE_PLL != 0  
  return *(partitions.partitionData[model]); 
#else 
  return tr.partitionData[model] ; 
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
  initReversibleGTR(&getTrHandle(), &(getPartitionsHandle()) , model); 
#else 
  initReversibleGTR(&getTrHandle(), model); 
#endif
}


void TreeAln::setFrequencies(const std::vector<double> &values, int model)
{
  assert( BoundsChecker::checkFrequencies(values) ) ;    
  auto& partition = getPartition(model); 
  memcpy( partition.frequencies, &(values[0]), partition.states * sizeof(double)); 
  initRevMat(model);   
}


void TreeAln::setRevMat(const std::vector<double> &values, int model)
{
  bool valuesOkay = BoundsChecker::checkRevmat(values); 
  if(not valuesOkay)
    {
      tout << "Problem with substitution parameter: "  << MAX_SCI_PRECISION << values << std::endl;  
      assert( valuesOkay ); 
    }

  auto& partition = getPartition(model) ; 
  memcpy(partition.substRates, &(values[0]), numStateToNumInTriangleMatrix(  partition.states) * sizeof(double)); 
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
  auto& p = getPartition(model); 
  p.alpha = alpha; 
  discretizeGamma(model); 
}


/**
   @brief makes the discrete categories for the gamma
   distribution. Has to be called, if alpha was updated.   

*/ 
void TreeAln::discretizeGamma(int model)
{
  auto& partition =  getPartition(model); 
  makeGammaCats(partition.alpha, partition.gammaRates, 4, tr.useMedian);
}


std::ostream& operator<< (std::ostream& out,  const TreeAln&  traln)
{
  auto tp = TreePrinter(true, false, false); 
  return out << tp.printTree(traln); 
}


std::vector<double> TreeAln::getRevMat(int model) const 
{
  auto result = std::vector<double>{} ; 
  auto& partition = getPartition(model); 
  double sum = 0; 
  for(nat i =0 ; i < numStateToNumInTriangleMatrix(partition.states); ++i)
    {
      result.push_back(partition.substRates[i]);
      sum += result[i]; 
    }

  for_each(result.begin(), result.end(), [&](double &v) { v /= sum ; }) ; 

  return result; 
}


std::vector<double> TreeAln::getFrequencies(int model) const
{
  std::vector<double> result; 
  auto& partition = getPartition(model) ; 
  for(int i = 0; i < partition.states; ++i) 
    result.push_back(partition.frequencies[i]); 
  return result; 
}



bool TreeAln::revMatIsImmutable(int model) const
{
#ifdef UNSURE
  assert(0); 
#endif
    
  auto& partition = getPartition(model); 
    
  return partition.states == 20 && partition.protModels != GTR; 
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

  return  getTrHandle().nodep[elem] ; 
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
  auto sum = double{0.}; 

#if UNSURE
  // TDOO is this really correct?
  assert(0);
#endif


#if HAVE_PLL == 0
  for(auto &p :partitions)
    {
      result += tr.fracchanges[p] * tr.partitionContributions[p];   
      sum += tr.partitionContributions[p]; 
    }
#else 
  for(auto &p : partitions)
    {
      auto& partition =  getPartition(p);
      result += partition.fracchange * partition.partitionContribution;  
      sum += partition.partitionContribution; 
    }
#endif

  return result / sum ; 
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
  nat length = partition.upper - partition.lower; 
#else 
  nat length = partition.width; 
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


void TreeAln::setNumberOfPartitions(nat numPart)
{
#if HAVE_PLL == 0
  tr.NumberOfModels = numPart; 
#else 
  partitions.numberOfPartitions = numPart; 
#endif
}


void TreeAln::setModelAssignment(int part, ProtModel model) 
{
  auto& pData = getPartition(part) ; 
  
  // tout << "setting model " << ProtModelFun::getName(model) << std::endl; 

  pData.protModels=int(model);
  initRevMat(part); 
}


ProtModel TreeAln::getModelAssignment(int part) const
{
  auto& pData = getPartition(part) ; 
  return ProtModel(pData.protModels); 
}
