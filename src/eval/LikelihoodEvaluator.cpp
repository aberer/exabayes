#include <algorithm>
#include <numeric>

#include "LikelihoodEvaluator.hpp"
#include "system/GlobalVariables.hpp"
#include "ArrayRestorer.hpp"
#include "model/Branch.hpp"
#include "comm/ParallelSetup.hpp"

#include "common.h"

using std::endl; 
using std::cout; 


using std::endl; 
using std::cout; 



LikelihoodEvaluator::LikelihoodEvaluator(const TreeAln &traln, ArrayPolicy* plcy , std::shared_ptr<ArrayReservoir> arrayReservoir, ParallelSetup* pl)
  : _arrayPolicy (plcy->clone()) 
  , _arrayOrientation(traln)
  , _arrayReservoir(arrayReservoir)
  , _plPtr(pl)
{
}


LikelihoodEvaluator::LikelihoodEvaluator( const LikelihoodEvaluator &rhs )
  : _debugTralnPtr(std::move(rhs._debugTralnPtr))
  , _verifyLnl(rhs._verifyLnl)
  , _arrayPolicy(std::move(rhs._arrayPolicy->clone()))
  , _arrayOrientation(rhs._arrayOrientation)
  , _arrayReservoir(rhs._arrayReservoir)
  , _plPtr(rhs._plPtr)
{
}


void swap(LikelihoodEvaluator &lhs, LikelihoodEvaluator  &rhs)
{
  using std::swap; 
  swap(lhs._debugTralnPtr, rhs._debugTralnPtr); 
  swap(lhs._arrayPolicy, rhs._arrayPolicy); 
  swap(lhs._arrayOrientation, rhs._arrayOrientation); 
  swap(lhs._arrayReservoir, rhs._arrayReservoir); // for consistency 
  swap(lhs._plPtr, rhs._plPtr); 
}


LikelihoodEvaluator& LikelihoodEvaluator::operator=(LikelihoodEvaluator rhs) 
{
  swap(*this, rhs); 
  return *this; 
}


void LikelihoodEvaluator::updatePartials (TreeAln &traln, nodeptr p)
{  
  pllInstance *tr = &traln.getTrHandle(); 
  partitionList *pr = &traln.getPartitionsHandle(); 

  if(isTip(p->number, tr->mxtips))
    return;

  /* the first entry of the traversal descriptor is always reserved for evaluate or branch length optimization calls,
     hence we start filling the array at the second entry with index one. This is not very nice and should be fixed 
     at some point */

  tr->td[0].count = 0;

  /* compute the traversal descriptor, which will include nodes-that-need-update descending the subtree  p */
  computeTraversal(tr, p, PLL_TRUE, pr->perGeneBranchLengths?pr->numberOfPartitions : 1);

  /* the traversal descriptor has been recomputed -> not sure if it really always changes, something to 
     optimize in the future */
  tr->td[0].traversalHasChanged = PLL_TRUE;

  /* if there is something to re-compute */

  if(tr->td[0].count > 0)
    {
      /* store execute mask in traversal descriptor */

      storeExecuteMaskInTraversalDescriptor(tr, pr);

      for(int i = 0; i < tr->td[0].count ; ++i)
	{
	  auto& tInfo = tr->td[0].ti[i];

	  for(nat model = 0 ; model < traln.getNumberOfPartitions(); ++model)
	    {
	      auto &partition = traln.getPartition(model).getHandle();

	      size_t            
		width  = partition.width;

	      if(tr->td[0].executeModel[model] && width > 0)
		{
		  size_t requiredLength = 0; 

		  // assert(not tr->saveMemory ) ; 

		  // determine how much memory we need 
		  if(tr->saveMemory)
		    {
		      size_t
			j,
			setBits = 0;                

		      auto x1_gap = &(partition.gapVector[tInfo.qNumber * partition.gapVectorLength]);
		      auto x2_gap = &(partition.gapVector[tInfo.rNumber * partition.gapVectorLength]);
		      auto x3_gap = &(partition.gapVector[tInfo.pNumber * partition.gapVectorLength]);
            
		      for(j = 0; j < (size_t)pr->partitionData[model]->gapVectorLength; j++)
			{              
			  x3_gap[j] = x1_gap[j] & x2_gap[j];
			  setBits += 
			    // bitcount_32_bit(x3_gap[j])
			    __builtin_popcount(x3_gap[j])
			    ; 
			
			  // std::cout << SyncOut() << "set bits is " << setBits << std::endl; 
			}

		      requiredLength = (width - setBits)  * Partition::maxCategories * partition.states * sizeof(double);            
		    }
		  else 
		    {
		      requiredLength  =  virtual_width( width ) * Partition::maxCategories * partition.states * sizeof(double);
		    }

		  auto x3_start = partition.xVector[tInfo.pNumber - tr->mxtips - 1]; 
		  auto availableLength = partition.xSpaceVector[(tInfo.pNumber - tr->mxtips - 1)]; 

		  // tout << "required = " << requiredLength << "\tavailableLength=" << availableLength << std::endl; 

		  /* new exabayes behaviour  */
		  if(requiredLength != availableLength)
		    {
		      int p_slot =  tInfo.pNumber - tr->mxtips - 1; 

		      if( x3_start)
			_arrayReservoir->deallocate(x3_start);

		      x3_start = _arrayReservoir->allocate( requiredLength);
		  
		      partition.xVector[p_slot] = x3_start; 
		      partition.xSpaceVector[p_slot] = requiredLength; 
		    }

		}
	    }
	}
      
      pllNewviewIterative(tr, pr, 0);
      
    }

  tr->td[0].traversalHasChanged = PLL_FALSE;
}


void LikelihoodEvaluator::evaluateLikelihood (TreeAln &traln , BranchPlain branch, bool fullTraversal, bool isDebugEval)
{
  pllInstance* tr = &(traln.getTrHandle()); 
  partitionList *pr = &(traln.getPartitionsHandle());
  nodeptr p = branch.findNodePtr(traln) ; 

  auto result = log_double::fromAbs(1.);
  
  nodeptr 
    q = p->back; 

  int
    numBranches = pr->numberOfPartitions;

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;          

  for(int i = 0; i < numBranches; i++)
    tr->td[0].ti[0].qz[i] =  q->z[i];

  tr->td[0].count = 1;

  if(fullTraversal)
    {
      assert(isTip(q->back->number, tr->mxtips));
      computeTraversal(tr, q, PLL_FALSE, numBranches);
    }
  else
    {
      if( needsRecomp(tr->useRecom, tr->rvec, p, tr->mxtips) )
	computeTraversal(tr, p, PLL_TRUE, numBranches);

      if(  needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips) )
	computeTraversal(tr, q, PLL_TRUE, numBranches);
    }

  // this function mostly exists for debug purposes. If it is not
  // called in debug mode, then this function really should only take
  // care of the evaluate part. All nodes up until the virtual root
  // should have been computed by some other function (preferably
  // evaluateSubtrees) already.
  if(not isDebugEval)
    {
      if(tr->td[0].count != 1)
      tout << "danger " << SHOW(tr->td[0].count) << " but expected 1" << std::endl; 
    }
  // assert( tr->td[0].count == 1 ); 

  storeExecuteMaskInTraversalDescriptor(tr, pr);

  tr->td[0].traversalHasChanged = PLL_TRUE;

  pllEvaluateIterative(tr, pr, PLL_FALSE); //PLL_TRUE
  
  auto lnls = traln.getPartitionLnls(); 

  // meh 
  {
    auto lnlRaw = std::vector<double>(lnls.size());
    for(auto i = 0u; i < lnls.size(); ++i)
      lnlRaw[i] = lnls[i].getRawLog();
  
    lnlRaw = _plPtr->getChainComm().allReduce(lnlRaw); 
  
    for(auto i = 0u; i < lnlRaw.size(); ++i)
      lnls[i] = log_double::fromLog(lnlRaw[i]);
  }

  traln.setPartitionLnls(lnls); 

  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    result *= lnls[i]; 

  tr->likelihood = result.getRawLog();    
  tr->td[0].traversalHasChanged = PLL_FALSE;
}


bool LikelihoodEvaluator::applyDirtynessToSubtree(TreeAln &traln, nat partId, const BranchPlain &branch) 
{
  if(traln.isTipNode(branch.getPrimNode()))
    return true; 

  nat id =  ArrayOrientation::nodeNum2ArrayNum( branch.getPrimNode() , traln.getNumberOfTaxa( ) ); 
  auto p = branch.findNodePtr(traln); 

  auto& partition = traln.getPartition(partId).getHandle();

// #if 1 
  bool isClean = partition.xSpaceVector[id] != 0  && _arrayOrientation.isCorrectNew(traln, partId, p ); 
// #else 
//   bool isClean = partition.xSpaceVector[id] != 0  && _arrayOrientation.isCorrect(partId, id, p->back->number); 
// #endif

  if( isClean)
    {
      p->x = 1; 
      p->next->x = 0 ; 
      p->next->next->x = 0; 
    }
  else 
    { 
      _arrayOrientation.setOrientation(partId, id, p->back->number); 

      assert(not traln.isTipNode(p->number)); 
      p->x = 0; 
      p->next->x = 1 ; 
      p->next->next->x = 0; 

      auto descs = traln.getDescendents(branch); 
      applyDirtynessToSubtree(traln, partId, descs.first.getInverted()); 
      applyDirtynessToSubtree(traln, partId, descs.second.getInverted());     
    }

  return isClean; 
}


void LikelihoodEvaluator::evaluate( TreeAln &traln, const BranchPlain &root, bool fullTraversal)    
{
  // tout << "evaluate with root " << root << " and " << ( fullTraversal ? "full" : "partial" ) << std::endl; 
  auto partitions = std::vector<nat>{}; 
  nat numPart = traln.getNumberOfPartitions(); 
  partitions.reserve(numPart);
  for(nat i = 0 ; i < numPart; ++i)
    partitions.push_back(i); 
  evaluatePartitionsWithRoot(traln, root, partitions, fullTraversal); 
}





void LikelihoodEvaluator::evaluateSubtrees(TreeAln& traln, const BranchPlain &root , const std::vector<nat> & partitions, bool fullTraversal)
{
  auto toExecute = std::vector<bool>(traln.getNumberOfPartitions(), false);   
  for(auto m : partitions)
    {
      auto &partition = traln.getPartition(m);
      toExecute[m] = partition.getWidth() > 0 ; 
    }
  traln.setExecModel(toExecute); 

  if(fullTraversal)
    {
      for(auto &elem : partitions)
	{
	  auto &partition = traln.getPartition(elem);
	  if(partition.getWidth() > 0 )
	    markPartitionDirty(traln, elem);
	}
    }

  for(auto elem : partitions )
    {
      auto &partition = traln.getPartition(elem);
      if(partition.getWidth() > 0 )
	{
	  evalSubtree(traln, elem, root, fullTraversal); 
	  evalSubtree(traln, elem, root.getInverted(), fullTraversal); 
	}
    }
}




void LikelihoodEvaluator::evaluatePartitionsDry(TreeAln &traln, const BranchPlain &root , const std::vector<nat>& partitions)
{
  // horrible hack!  

  auto widthBackup = std::vector<nat>{}; 
  for(auto &p : partitions)
    {
      auto &ph = traln.getPartition(p).getHandle();
      widthBackup.push_back(ph.width);
      ph.width = 0; 
    }

  evaluatePartitionsWithRoot(traln, root, partitions, false);
  
  for(nat i = 0; i < widthBackup.size(); ++i)
    {
      auto p = partitions[i]; 
      auto &ph = traln.getPartition(p).getHandle();
      ph.width = widthBackup[i]; 
    }
}




void LikelihoodEvaluator::evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain &root , const std::vector<nat>& partitions, bool fullTraversal)    
{
  // tout << "evaluatePartitions with root " << root << " and " << ( fullTraversal ? "full" : "partial" )  << " and partitons " << partitions << std::endl; 
  nat numPart = traln.getNumberOfPartitions();
  auto perPartitionLH = traln.getPartitionLnls();

  evaluateSubtrees(traln, root, partitions, fullTraversal);

  auto execute = std::vector<bool>(traln.getNumberOfPartitions(), false); 
  for (auto elem : partitions)
    {
      auto &partition = traln.getPartition(elem); 
      execute[elem] = partition.getWidth() > 0; 
    }
  traln.setExecModel(execute); 

  bool changedOrientation  = std::any_of(begin(execute), end(execute), [](bool elem){return elem ; });

  auto start = root.findNodePtr(traln);
  if( changedOrientation && not ( ( start->x == 1 || traln.isTipNode(start)   )
	   && ( start->back->x == 1 || traln.isTipNode(start->back) )  ) )
    {
      tout << "at evaluation at " << root << ": " ; 

      if( not ( start->x == 1 || traln.isTipNode(start)  ))
	{
	  tout << " problem with " << start->number << std::endl; 
	  assert(0); 
	}
    }
  evaluateLikelihood(traln, root, false, false); 

  auto pLnl = traln.getPartitionLnls();
  for(auto m : partitions )
    perPartitionLH[m] = pLnl[m]; 
  traln.setPartitionLnls(perPartitionLH); 

  traln.setLikelihood ( std::accumulate(begin(perPartitionLH), end(perPartitionLH), log_double::fromAbs(1.0), std::multiplies<log_double>()) ) ; 
  traln.setExecModel(std::vector<bool>(numPart, true));

#ifdef DEBUG_LNL_VERIFY
  expensiveVerify(traln, root ,traln.getLikelihood());   
#endif
}


/** 
    NOTE: the fullTraversal parameter is only there because of the searching (HACK)
 */ 
void LikelihoodEvaluator::evalSubtree(TreeAln  &traln, nat partId, const BranchPlain &evalBranch, bool fullTraversal)   
{
  // TODO do we have a cleaner partition for this break? It is
  // necessary, s.t. in the -Q case, we do not do all that traverals
  // (prepareForEvaluation) for partitions that do not have data anyway
  auto &partition = traln.getPartition(partId);
  if(partition.getWidth() == 0 )
    return; 
  
  // tout << "eval subtree for " << partId << std::endl; 

  if(traln.isTipNode(evalBranch.getPrimNode()))
    {
      // tout << evalBranch << " is a tip node " << std::endl; 
      return; 
    }

  _arrayPolicy->prepareForEvaluation(traln, evalBranch, partId, _arrayOrientation, *_arrayReservoir); 
  applyDirtynessToSubtree(traln, partId, evalBranch);
  
  if( 
#ifndef LNL_PRINT_DEBUG
     false && 
#endif
     true  
      )
    {
      if(evalBranch.findNodePtr(traln)->x != 1)
	{
	  tout << "computing "; 
	  debugPrintToComputeHelper(traln, evalBranch); 
	  tout << std::endl; 
	}
      else 
	{
	  tout << "nothing to compute" << std::endl; 
	}      
    }
  
  // abort, if everything is okay 
  if(evalBranch.findNodePtr(traln)->x == 1)
    return; 

  auto execute = std::vector<bool>(traln.getNumberOfPartitions(), false); 
  execute[partId] = true; 
  traln.setExecModel(execute); 

  updatePartials(traln, evalBranch.findNodePtr(traln)	); 

  // assert stuff 
  auto desc = traln.getDescendents(evalBranch); 
  auto p = desc.first.getInverted().findNodePtr(traln); 
  assert(p->x || traln.isTipNode(p->number)); 
  p = desc.second.getInverted().findNodePtr(traln); 
  assert(p->x || traln.isTipNode(p->number)); 

  p =  evalBranch.findNodePtr(traln) ; 

  p->x = 1; 
  p->next->x  = 0; 
  p->next->next->x = 0; 
}


void LikelihoodEvaluator::disorientDebugHelper(TreeAln &traln, const BranchPlain& root  )
{
  auto p = root.findNodePtr(traln); 

  if( traln.isTipNode(p))
    return ; 

  p->x = 0; 
  p->next->x = 1; 
  p->next->next->x = 0; 

  auto descs = traln.getDescendents(root); 
  disorientDebugHelper(traln, descs.first.getInverted()); 
  disorientDebugHelper(traln, descs.second.getInverted()); 
}


void LikelihoodEvaluator::disorientDebug(TreeAln &traln, const BranchPlain& root  )
{
  disorientDebugHelper(traln,root); 
  disorientDebugHelper(traln, root.getInverted()); 
}


void LikelihoodEvaluator::expensiveVerify(TreeAln &traln, BranchPlain root, log_double toVerify )
{
  // std::cout << "orig-lnl=" << traln.getTrHandle().likelihood << std::endl;  
  debugPrint = 0 ; 
  
  _debugTralnPtr->clearMemory();
  *_debugTralnPtr = traln; 

  disorientDebug(*_debugTralnPtr, root); 
  evaluateLikelihood(*_debugTralnPtr, root, false, true);   

  auto verifiedLnl =  _debugTralnPtr->getLikelihood(); 

  if( fabs(log_double(verifiedLnl / toVerify ).toAbs())  - 1.0  > ACCEPTED_LIKELIHOOD_EPS)
    {
      for(nat j = 0; j < traln.getNumberOfPartitions() ; ++j )
	{
	  auto& partition = _debugTralnPtr->getPartition(j).getHandle(); 
	  auto sc = partition.globalScaler; 
	  tout << "scCorr=" ; 
	  nat ctr = 0; 
	  for(nat i = 0; i < 2 * traln.getNumberOfTaxa(); ++i)
	    {
	      if(sc[i] != 0 )
		tout << ctr <<  "=" << sc[i] << ","; 
	      ++ctr; 
	    }
	  tout << std::endl; 

	  auto& partition2 = traln.getPartition(j).getHandle(); 
	  sc = partition2.globalScaler; 
	  tout << "scReal=" ; 
	  ctr = 0; 
	  for(nat i = 0; i < 2 * traln.getNumberOfTaxa(); ++i)
	    {
	      if(sc[i] != 0)
		{
		  tout << ctr << "=" << sc[i] << " "; 
		}
	      ++ctr; 
	    }
	  tout << std::endl; 
	} 

      tout << "WARNING: found in expensive evaluation: likelihood difference is " 
	   << " TODO "
	   // <<  std::setprecision(8) <<   fabs (verifiedLnl - toVerify )
	   << " (with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")\n"
	   << "partitionLnls= " << traln.getPartitionLnls() << "\n"
	   << "verifiedPartLnls="  << _debugTralnPtr->getPartitionLnls() << "\n"; 

      evaluateLikelihood(traln, traln.getAnyBranch(), true, true); 
      tout  << "after full traversal on orig=" << traln.getTrHandle().likelihood << std::endl; 
      assert(0); 
    }  
  // assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS); 
  
  debugPrint = 1;  

  // std::cout << SyncOut() << "VERIFIED " << MAX_SCI_PRECISION << verifiedLnl << SOME_FIXED_PRECISION << std::endl; 

  _debugTralnPtr->clearMemory();
}


#ifdef DEBUG_LNL_VERIFY
void LikelihoodEvaluator::setDebugTraln(std::shared_ptr<TreeAln> debugTraln)
{
  _verifyLnl = true; 
  _debugTralnPtr = debugTraln; 
}
#endif


void LikelihoodEvaluator::markDirty(const TreeAln &traln, nat nodeId)
{
  for(nat i = 0; i < traln.getNumberOfPartitions() ;++i)
    _arrayOrientation.setOrientation(i, ArrayOrientation::nodeNum2ArrayNum(nodeId,traln.getNumberOfTaxa()), INVALID); 
} 


void LikelihoodEvaluator::debugPrintToComputeHelper(const TreeAln &traln, const BranchPlain &root )
{
  if(traln.isTipNode(root.getPrimNode()))
    return; 

  auto desc = traln.getDescendents(root); 
  
  if( root.findNodePtr(traln)->x == 0)
    {
      tout << root.getPrimNode() <<  "=(" << desc.first.getSecNode( )<<  "+"  << desc.second.getSecNode()  << "),"  ; 
      debugPrintToComputeHelper(traln, desc.first.getInverted()); 
      debugPrintToComputeHelper(traln, desc.second.getInverted()); 
    }
}


// prints what will be computed 
void LikelihoodEvaluator::debugPrintToCompute(const TreeAln &traln, const BranchPlain &root)
{
  tout << "to compute: 0=(" <<   root.getPrimNode() << "+" << root.getSecNode() << ")," ; 
  debugPrintToComputeHelper(traln, root); 
  debugPrintToComputeHelper(traln, root.getInverted()); 
}


void LikelihoodEvaluator::markDirty(const TreeAln &traln, nat partitionId, nat nodeId) 
{
  nat id = nodeId - traln.getNumberOfTaxa()  - 1; 
  _arrayOrientation.setOrientation(partitionId,id, INVALID); 
} 

void LikelihoodEvaluator::markPartitionDirty(const TreeAln &traln, nat partition )
{
  _arrayOrientation.setPartitionInvalid(partition);
}

void LikelihoodEvaluator::accountForRejection(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat> &invalidNodes)
{
  _arrayPolicy->accountForRejection(traln,partitions, invalidNodes, _arrayOrientation,*_arrayReservoir); 
}

void LikelihoodEvaluator::freeMemory()
{ 
  _arrayPolicy->freeMemory(*_arrayReservoir) ; 
}  

void LikelihoodEvaluator::imprint(const TreeAln &traln) 
{ 
  _arrayPolicy->imprint(traln, _arrayOrientation);  
}



void  LikelihoodEvaluator::evaluateSubtrees(TreeAln& traln, const BranchPlain &root, bool fullTraversal)
{
  auto parts = std::vector<nat>();
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    parts.push_back(i);
  evaluateSubtrees(traln, root, parts, fullTraversal); 
} 
