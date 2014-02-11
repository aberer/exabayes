#include <sstream>
#include <cassert>
#include <cstring>
#include <unordered_set>

#include "RateHelper.hpp"
#include "tree-init/TreeResource.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/AbstractParameter.hpp"
#include "TreeRandomizer.hpp"
#include "axml.h"
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"
#include "BoundsChecker.hpp"
#include "TreePrinter.hpp"
#include "tree-init/TreeInitializer.hpp"


TreeAln::TreeAln(nat numTax)
  : _mode(RunModes::NOTHING)
  , _hasLnlArrays(false)
  , _hasAlignment{false}
  , _hasTopology{false}
{
  memset(&_tr,0,sizeof(tree)); 
  createStandardTree(numTax);
}


TreeAln::TreeAln( const TreeAln& rhs)
  : _mode(rhs._mode)
  , _taxa(rhs._taxa)
  , _hasLnlArrays(false)
  , _hasAlignment(false)
  , _hasTopology(false)
{
  memset(&_tr,0,sizeof(tree));
  createStandardTree(rhs.getNumberOfTaxa());

  if(rhs._hasAlignment)
    {
      // tout << "copying alignment " << std::endl; 
      auto &&ti = TreeInitializer(std::unique_ptr<InitializationResource>(new TreeResource(&rhs)));
      ti.initializeWithAlignmentInfo(*this, rhs._mode);
      _hasAlignment = true; 
      copyAlnModel(rhs); 
    }

  // tout << "copying topology" << std::endl; 
  if(rhs._hasTopology)
    {
      copyTopologyAndBl(rhs);
      _hasTopology = true; 
    }
}


// this was problematic with gcc-4.6
// TreeAln::TreeAln(TreeAln &&rhs)
// {
//   // memset(&_tr, 0, sizeof(tree)); 
//   // createStandardTree(rhs.getNumberOfTaxa()); 

//   // assert(0); 
//   // tout << "MOVE tree " << &rhs << std::endl; 
//   // std::swap(*this, rhs); 
//   // using std::swap; 
//   swap(*this, rhs);
// }


TreeAln& TreeAln::operator=(TreeAln rhs)
{
  // assert(0); 
  // tout << "using OPERATOR= on " << this <<  " to  copy rhs=" << &rhs  << std::endl; 
  // using std::swap; 
  swap(*this, rhs);
  return *this; 
}



void swap(TreeAln& lhs, TreeAln& rhs )
{
  using std::swap; 
#if HAVE_PLL != 0
  swap(lhs._partitions, rhs._partitions); 
#endif 
  swap(lhs._taxa, rhs._taxa); 
  swap(lhs._mode, rhs._mode); 
  swap(lhs._tr, rhs._tr); 
  
  swap(lhs._hasLnlArrays,  rhs._hasLnlArrays); 
  swap(lhs._hasAlignment ,rhs._hasAlignment );
  swap(lhs._hasTopology , rhs._hasTopology); 
}




TreeAln::~TreeAln()
{  
  // this is slightly crazy, but we'll have to do the cleanup
  // somewhere.

  if(_hasAlignment)
    {
      nat numTax = getNumberOfTaxa(); 

      delete [] _tr.zqr; 
      delete [] _tr.currentZQR; 
      delete [] _tr.currentLZR; 
      delete [] _tr.currentLZQ; 
      delete [] _tr.currentLZS; 
      delete [] _tr.currentLZI; 
      delete [] _tr.lzs; 
      delete [] _tr.lzq; 
      delete [] _tr.lzr; 
      delete [] _tr.lzi; 
      delete [] _tr.coreLZ; 
      delete [] _tr.curvatOK; 
      
      delete [] _tr.partitionSmoothed ; 
      delete [] _tr.partitionConverged; 

      for(nat i = 0; i < numTax ; ++i)
	{
	  delete [] _tr.td[0].ti[i].qz; 
	  delete [] _tr.td[0].ti[i].rz; 
	}

      for(nat i = 0; i < numTax + 3 * (numTax - 1)  ; ++i )
	{
	  auto node = _tr.nodeBaseAddress + i; 
	  delete [] node->z; 
	}
    }

  if(_tr.aliaswgt != NULL)
    exa_free(_tr.aliaswgt); 
#if HAVE_PLL == 0
  if(_tr.perPartitionLH != NULL)
    exa_free(_tr.perPartitionLH); 
#endif
  if(_tr.rateCategory != NULL)
    exa_free(_tr.rateCategory); 
  if(_tr.patrat != NULL)
    exa_free(_tr.patrat);
  if(_tr.patratStored != NULL)
    exa_free(_tr.patratStored);
  if(_tr.lhs != NULL)
    exa_free(_tr.lhs);  
  if(_tr.yVector != NULL)
    exa_free(_tr.yVector); 

  if(_tr.tree_string != NULL)
    exa_free(_tr.tree_string);
  if(_tr.tree0 != NULL)
    exa_free(_tr.tree0);
  if(_tr.tree1 != NULL)
    exa_free(_tr.tree1);
  
  if(_tr.constraintVector != NULL)
    exa_free(_tr.constraintVector); 
  if(_tr.nodeBaseAddress != NULL)
    exa_free(_tr.nodeBaseAddress);

  // free parsimony related stuff
  exa_free(_tr.parsimonyScore);
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
      auto& partition = getPartition(i); 
      if(partition.parsVect != NULL)
	exa_free(partition.parsVect);
    }
  exa_free(_tr.ti); 

  exa_free(_tr.td[0].ti);    
  exa_free(_tr.td[0].executeModel);    
  exa_free(_tr.td[0].parameterValues);    

  exa_free(_tr.nodep); 

  int numPart = getNumberOfPartitions();
  
  // tout << "have "  << numPart << std::endl; 

  if(_hasAlignment)
    {
      for(int i = 0; i < numPart ;++i)
	{
	  auto& partition = getPartition(i); 
	  
	  if( _tr.saveMemory)
	    {
	      if(partition.gapColumn != 0 )
		delete [] partition.gapColumn; 
	      if(partition.gapVector != 0)
		delete [] partition.gapVector; 
	    }

#if HAVE_PLL ==  0
	  exa_free(partition.symmetryVector);
	  exa_free(partition.frequencyGrouping); 
#endif

	  exa_free(partition.xVector); 

	  exa_free(partition.frequencies);
	  exa_free(partition.empiricalFrequencies); 
	  exa_free(partition.tipVector); 
	  exa_free(partition.perSiteRates); 
	  exa_free(partition.gammaRates); 
      
	  for(int j = 1; j < _tr.mxtips + 1 ; ++j)
	    {
	      if(partition.yVector[j] != NULL)
		exa_free(partition.yVector[j]);
	    }
	  exa_free(partition.yVector);
	  exa_free(partition.xSpaceVector); 

#if HAVE_PLL == 0
	  exa_free(partition.sumBuffer); 
#endif

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


      if(getNumberOfPartitions() > 0)
	{
#if HAVE_PLL != 0 
	  for(nat i = 0; i < getNumberOfPartitions(); ++i)
	    {
	      exa_free(_partitions.partitionData[i]);
	    }
	  exa_free(_partitions.partitionData); 
#else 
	  exa_free(_tr.partitionData);
#endif

	}

    }

#if HAVE_PLL == 0  
  exa_free(_tr.executeModel) ;  
  exa_free(_tr.partitionContributions); 
  exa_free(_tr.fracchanges); 
#endif
}


void TreeAln::createStandardTree(nat numTax)
{
  _tr.likelihood =  0 ; 
  _tr.doCutoff = TRUE;
  _tr.secondaryStructureModel = SEC_16; /* default setting */
  _tr.searchConvergenceCriterion = FALSE;
  _tr.rateHetModel = GAMMA; 
  _tr.multiStateModel  = GTR_MULTI_STATE;
#if HAVE_PLL == 0 
  _tr.useGappedImplementation = FALSE;
  _tr.saveBestTrees          = 0;
  _tr.numBranches = getNumberOfPartitions();
#else 
  _partitions.perGeneBranchLengths = TRUE; 
#endif

  _tr.manyPartitions = FALSE;
  _tr.saveMemory = FALSE;
  _tr.categories             = 25;
  _tr.grouped = FALSE;
  _tr.constrained = FALSE;
  _tr.gapyness               = 0.0; 
  _tr.useMedian = FALSE;
  _tr.mxtips = numTax; 

  setNumberOfPartitions(0);

  TreeInitializer::setupTheTree(_tr);
}


void TreeAln::clearMemory(ArrayReservoir &arrayReservoir)
{
  for(nat i = 0; i < getNumberOfPartitions() ; ++i)
    {
      auto& partition = getPartition(i);
      for(nat j = 0; j < getNumberOfTaxa() ; ++j)
	{
	  if(partition.xSpaceVector[j] != 0 )
	    {
	      if(partition.xVector[j])
		{
		  arrayReservoir.deallocate(partition.xVector[j]); 
		  partition.xVector[j] = nullptr; 
		  partition.xSpaceVector[j] = 0; 
		}
	    }
	}
    }
}


void TreeAln::clipNodeDefault(nodeptr p, nodeptr q )
{
  clipNode(p,q);  
  // use a problematic branch length to ensure, that it was overridden
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    p->z[i] = p->back->z[i] = 0.9; 
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
  for(int i = 1; i < _tr.mxtips+1; ++i)
    {
      nodeptr p = _tr.nodep[i]; 
      p->back = NULL; 
    }

  for(int i = _tr.mxtips + 1; i < 2 * _tr.mxtips; ++i)
    unlinkNode(_tr.nodep[i]); 
}


nodeptr TreeAln::getUnhookedNode(int number)
{  
  nodeptr p = _tr.nodep[number]; 

  if(isTip(number, _tr.mxtips) )
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


void TreeAln::copyTopologyAndBl(const TreeAln &rhs) 
{
  assert(&rhs != this); 

  this->unlinkTree();
  int mxtips = rhs.getTrHandle().mxtips ;
  auto &rhsTree = rhs.getTrHandle(); 
  auto &thisTree = getTrHandle();

  // if this works, it is one of the most hackney things, I've ever done... 
  for(int i = mxtips+1 ; i < 2* _tr.mxtips-1; ++i)
    {
      nodeptr rhsNode = rhsTree.nodep[i],
	lhsNode = thisTree.nodep[i];       
      hookup(lhsNode, lhsNode -   (rhsNode - rhsNode->back), rhsNode->z, rhs.getNumberOfPartitions() ) ;             

      lhsNode = lhsNode->next; 
      rhsNode = rhsNode->next; 
      hookup(lhsNode, lhsNode -   (rhsNode - rhsNode->back), rhsNode->z, rhs.getNumberOfPartitions()) ;       

      lhsNode = lhsNode->next; 
      rhsNode = rhsNode->next; 
      hookup(lhsNode, lhsNode  -   (rhsNode - rhsNode->back), rhsNode->z, rhs.getNumberOfPartitions()) ; 
    }
}


void TreeAln::copyAlnModel(const TreeAln& rhs)
{
  assert(&rhs != this); 

 // copy partition parameters 
  for(nat i = 0; i < rhs.getNumberOfPartitions(); ++i)
    {
      auto& partitionRhs = rhs.getPartition(i); 
      auto& partitionLhs = this->getPartition(i); 

      for(int j = 0; j < partitionRhs.states; ++j )
	partitionLhs.frequencies[j] = partitionRhs.frequencies[j]; 
      
      for(nat j = 0 ; j < RateHelper::numStateToNumInTriangleMatrix(partitionRhs.states); ++j)
	partitionLhs.substRates[j] = partitionRhs.substRates[j]; 

      partitionLhs.protModels = partitionRhs.protModels; 	
      partitionLhs.alpha = partitionRhs.alpha; 

      partitionLhs.protFreqs = partitionRhs.protFreqs; 

      initRevMat(i);
      discretizeGamma(i);       
    }
}


nat TreeAln::getNumberOfPartitions() const
{
#if HAVE_PLL != 0
  return _partitions.numberOfPartitions; 
#else 
  return _tr.NumberOfModels; 
#endif
}


pInfo& TreeAln::getPartition(nat model)  const
{
  assert(model < getNumberOfPartitions()); 
  
#if HAVE_PLL != 0  
  return *(_partitions.partitionData[model]); 
#else 
  return _tr.partitionData[model] ; 
#endif
}


std::vector<bool> TreeAln::getExecModel() const 
{ 
  std::vector<bool> result; 
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0  
      result.push_back(_partitions.partitionData[i]->executeModel); 
#else  
      result.push_back(_tr.executeModel[i]); 
#endif
    }
  return result; 
}
 
void TreeAln::setExecModel(const std::vector<bool>  &modelInfo)
{
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0 
      _partitions.partitionData[i]->executeModel =  modelInfo[i] ?  true : false; 
#else 
      _tr.executeModel[i] = modelInfo[i]; 
#endif
    }
}

std::vector<double> TreeAln::getPartitionLnls() const
{
  std::vector<double> result; 
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0 
      result.push_back(_partitions.partitionData[i]->partitionLH); 
#else 
      result.push_back(_tr.perPartitionLH[i]); 
#endif
    }
  return result; 
}
 
void TreeAln::setPartitionLnls(const std::vector<double> partitionLnls)  
{
  for(nat i = 0; i < getNumberOfPartitions(); ++i)
    {
#if HAVE_PLL != 0
      _partitions.partitionData[i]->partitionLH = partitionLnls[i]; 
#else 
      _tr.perPartitionLH[i]  = partitionLnls[i]; 
#endif
    }
} 


void TreeAln::initRevMat(nat model)
{
#if HAVE_PLL != 0
  initReversibleGTR(&getTrHandle(), &(getPartitionsHandle()) , model); 
#else 
  initReversibleGTR(&getTrHandle(), model); 
#endif
}

void TreeAln::setRevMat(const std::vector<double> &values, nat model )
{
  bool valuesOkay = BoundsChecker::checkRevmat(values); 
  if(not valuesOkay)
    {
      tout << "Problem with substitution parameter: "  << MAX_SCI_PRECISION << values << std::endl;  
      assert( valuesOkay ); 
    }

  auto& partition = getPartition(model) ; 

  assert(partition.dataType != AA_DATA || partition.protModels == GTR); 

  memcpy(partition.substRates, values.data(), RateHelper::numStateToNumInTriangleMatrix(  partition.states) * sizeof(double)); 
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


void TreeAln::setAlpha(double alpha,  nat model)
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
void TreeAln::discretizeGamma(nat model)
{
  auto& partition =  getPartition(model); 
  makeGammaCats(partition.alpha, partition.gammaRates, 4, _tr.useMedian);
}


std::ostream& operator<< (std::ostream& out,  const TreeAln&  traln)
{
  auto tp = TreePrinter(true, false, false); 
  return out << tp.printTree(traln); 
}


std::vector<double> TreeAln::getRevMat(nat model) const 
{
  auto result = std::vector<double>{} ; 
  auto& partition = getPartition(model); 
  assert(partition.dataType != AA_DATA || partition.protModels == GTR); 

  for(nat i = 0; i < RateHelper::numStateToNumInTriangleMatrix(partition.states); ++i)
    result.push_back(partition.substRates[i]);
  
  RateHelper::convertToSum1(result);

  return result; 
}


void TreeAln::setFrequencies(const std::vector<double> &values, nat model)
{
  assert( BoundsChecker::checkFrequencies(values) ) ;    
  auto& partition = getPartition(model); 
  assert( partition.dataType != AA_DATA || partition.protFreqs == TRUE ); 

  memcpy( partition.frequencies, values.data(), partition.states * sizeof(double)); 
  initRevMat(model);   
}


std::vector<double> TreeAln::getFrequencies(nat model) const
{
  auto result = std::vector<double> {}; 
  auto& partition = getPartition(model) ; 

  assert( partition.dataType != AA_DATA || partition.protFreqs == TRUE ); 
  
  for(int i = 0; i < partition.states; ++i) 
    result.push_back(partition.frequencies[i]); 
  return result; 
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

  for(auto &p :partitions)
    {
#if HAVE_PLL == 0
      result += _tr.fracchanges[p] * _tr.partitionContributions[p];   
      assert(_tr.partitionContributions[p] > 0 && _tr.fracchanges[p] > 0  ); 
      sum += _tr.partitionContributions[p]; 
#else 
      auto& partition =  getPartition(p);
      result += partition.fracchange * partition.partitionContribution;  
      sum += partition.partitionContribution; 
      assert(partition.partitionContribution > 0 && partition.fracchange > 0  ); 
#endif
    }

  return result / sum   ; 
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
  return BranchPlain(_tr.nodep[1]->number, _tr.nodep[1]->back->number); 
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
  _tr.NumberOfModels = numPart; 
#else 
  _partitions.numberOfPartitions = numPart; 
#endif
}


void TreeAln::setProteinModel(int part, ProtModel model) 
{
  auto& pData = getPartition(part) ; 
  pData.protModels=int(model);
  initRevMat(part); 
}


ProtModel TreeAln::getProteinModel(int part) const
{
  auto& pData = getPartition(part) ; 
  return ProtModel(pData.protModels); 
}


// HACK, don't use 
void TreeAln::setBranchUnchecked(const BranchLength &bl)
{
  auto p =  bl.findNodePtr(*this);
  p->z[0] = bl.getLength();
  p->back->z[0] = bl.getLength();
}


void TreeAln::clearMemory()
{
  for(nat i = 0; i < getNumberOfPartitions() ; ++i )
    {
      auto &partition = getPartition(i);
      for(nat j = 0; j < getNumberOfTaxa(); ++j)
	if(partition.xSpaceVector[j] != 0 )
	  {
	    exa_free(partition.xVector[j] );
	  }
    }
} 
