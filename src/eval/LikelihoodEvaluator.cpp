#include "LikelihoodEvaluator.hpp"
#include "GlobalVariables.hpp"
#include "ArrayRestorer.hpp"
#include "Branch.hpp"


LikelihoodEvaluator::LikelihoodEvaluator(const TreeAln &traln, ArrayPolicy* plcy )
  : arrayPolicy (plcy->clone()) 
  , arrayOrientation(traln)
{
}

LikelihoodEvaluator::LikelihoodEvaluator( const LikelihoodEvaluator &rhs )
  : debugTraln(std::move(rhs.debugTraln))
  , verifyLnl(rhs.verifyLnl)
  , arrayPolicy(std::move(rhs.arrayPolicy->clone()))
  , arrayOrientation(rhs.arrayOrientation)
{
}


LikelihoodEvaluator::LikelihoodEvaluator( LikelihoodEvaluator &&rhs )
  : debugTraln(std::move(rhs.debugTraln))
  , verifyLnl(rhs.verifyLnl)
  , arrayPolicy(std::move(rhs.arrayPolicy->clone()))
  , arrayOrientation(std::move(rhs.arrayOrientation))
{
}

void swap(LikelihoodEvaluator &lhs, LikelihoodEvaluator  &rhs)
{
  using std::swap; 
  swap(lhs.debugTraln, rhs.debugTraln); 
  swap(lhs.arrayPolicy, rhs.arrayPolicy); 
  swap(lhs.arrayOrientation, rhs.arrayOrientation); 
}

LikelihoodEvaluator& LikelihoodEvaluator::operator=(LikelihoodEvaluator rhs) 
{
  swap(*this, rhs); 
  return *this; 
}


bool LikelihoodEvaluator::applyDirtynessToSubtree(TreeAln &traln, nat partId, const BranchPlain &branch) 
{
  if(traln.isTipNode(branch.getPrimNode()))
    return true; 

  nat id = branch.getPrimNode() - traln.getNumberOfTaxa( ) -1 ; 

  auto p = branch.findNodePtr(traln); 

  auto& partition = traln.getPartition(partId);
  bool isClean = arrayOrientation.isCorrect(partId, id, p->back->number) && partition.xSpaceVector[id] != 0; 

  if( isClean)
    {
      p->x = 1; 
      p->next->x = 0 ; 
      p->next->next->x = 0; 
    }
  else 
    { 
      arrayOrientation.setOrientation(partId, id, p->back->number); 

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


void LikelihoodEvaluator::evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain &root , const std::vector<nat>& partitions, bool fullTraversal)    
{
  // tout << "evaluatePartitions with root " << root << " and " << ( fullTraversal ? "full" : "partial" )  << " and partitons " << partitions << std::endl; 

  nat numPart = traln.getNumberOfPartitions(); 

  auto perPartitionLH = traln.getPartitionLnls();

  auto toExecute = std::vector<bool>(numPart, false);   
  for(auto m : partitions)
    toExecute[m] = true; 
  traln.setExecModel(toExecute); 

  if(fullTraversal)
    {
      for(auto &elem : partitions)
	markPartitionDirty(traln, elem);
    }

  for(auto elem : partitions )
    {
      evalSubtree(traln, elem, root); 
      evalSubtree(traln, elem, root.getInverted()); 
    }

  auto execute = std::vector<bool>(traln.getNumberOfPartitions(), false); 
  for (auto elem : partitions)
    execute[elem] = true; 
  traln.setExecModel(execute); 

  exa_evaluateGeneric(traln, root);

  auto pLnl = traln.getPartitionLnls();
  for(auto m : partitions )
    perPartitionLH[m] = pLnl[m]; 
  traln.setPartitionLnls(perPartitionLH); 

  traln.getTrHandle().likelihood = std::accumulate(perPartitionLH.begin(), perPartitionLH.end(), 0.); 
  traln.setExecModel(std::vector<bool>(numPart, true));

#ifdef DEBUG_LNL_VERIFY
  expensiveVerify(traln, root ,traln.getTrHandle().likelihood);   
#endif
}


void LikelihoodEvaluator::evalSubtree(TreeAln  &traln, nat partition, const BranchPlain &evalBranch)   
{
  // tout << "eval subtree for " << partition << std::endl; 

  if(traln.isTipNode(evalBranch.getPrimNode()))
    {
      // tout << evalBranch << " is a tip node " << std::endl; 
      return; 
    }

  arrayPolicy->prepareForEvaluation(traln, evalBranch,  partition , arrayOrientation); 
  applyDirtynessToSubtree(traln, partition, evalBranch);
  
  if( 
     false && 
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
  execute[partition] = true; 
  traln.setExecModel(execute); 
  
  coreEvalSubTree(traln,evalBranch); 

  // assert stuff 
  auto desc = traln.getDescendents(evalBranch); 
  auto p = desc.first.getInverted().findNodePtr(traln); 
  assert(p->x || traln.isTipNode(p->number)); 
  p = desc.second.getInverted().findNodePtr(traln); 
  assert(p->x || traln.isTipNode(p->number)); 

  p =  evalBranch.findNodePtr(traln) ; 
  // tout << "setting root from " << p->number << " to "<< p->back->number << std::endl; 

  p->x = 1; 
  p->next->x  = 0; 
  p->next->next->x = 0; 
}


void LikelihoodEvaluator::exa_evaluateGeneric(TreeAln &traln, const BranchPlain& root )
{
  // tout << "calling exa_evaluate "  << std::endl; 

  auto start = root.findNodePtr(traln);

  if(not ( ( start->x == 1 || traln.isTipNode(start)   )
	   && ( start->back->x == 1 || traln.isTipNode(start->back) )  ) )
    {
      tout << "at evaluation at " << root << ": " ; 

      if( not ( start->x == 1 || traln.isTipNode(start)  ))
	{
	  tout << " problem with " << start->number << std::endl; 
	  assert(0); 
	}
    }

#if HAVE_PLL != 0
  evaluateGeneric(&traln.getTrHandle(), &traln.getPartitionsHandle(), start, FALSE); 
#else 
  evaluateGeneric(&traln.getTrHandle(), start, FALSE); 
#endif  
}


void LikelihoodEvaluator::coreEvalSubTree(TreeAln& traln, const BranchPlain &root)
{
  auto p = root.findNodePtr(traln); 
#if HAVE_PLL != 0
  newviewGeneric(&traln.getTrHandle(), &traln.getPartitionsHandle(), p, TRUE   ); 
#else 
  newviewGeneric(&traln.getTrHandle(), p, TRUE  ); 
#endif 
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



void LikelihoodEvaluator::expensiveVerify(TreeAln &traln, const BranchPlain& root, double toVerify)
{
  debugPrint = 0 ; 

  debugTraln->copyModel(traln); 

  disorientDebug(*debugTraln, root); 
  
#if HAVE_PLL != 0
  evaluateGeneric(&debugTraln->getTrHandle(), &debugTraln->getPartitionsHandle(), root.findNodePtr(*debugTraln), FALSE); 
#else 
  evaluateGeneric(&debugTraln->getTrHandle() , root.findNodePtr(*debugTraln), FALSE); 
#endif  

  double verifiedLnl =  debugTraln->getTrHandle().likelihood; 

  if(fabs (verifiedLnl - toVerify ) > ACCEPTED_LIKELIHOOD_EPS)
    {
#if 1 
      auto& partition = debugTraln->getPartition(0); 
      auto sc = partition.globalScaler; 
      tout << "scCorr=" ; 
      nat ctr = 0; 
      for(nat i = 0; i < 2 * traln.getNumberOfTaxa(); ++i)
	{
	  if(sc[i] != 0 )
	    {
	      tout << ctr <<  "=" << sc[i] << ","; 
	    }
	  ++ctr; 
	}
      tout << std::endl; 

      auto& partition2 = traln.getPartition(0); 
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
#endif

      std::cerr  << "WARNING: found in expensive evaluation: likelihood difference is " 
		 <<  std::setprecision(8) <<   fabs (verifiedLnl - toVerify )
		 << " (with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")" << std::endl; 

      tout << "partitionLnls= " << traln.getPartitionLnls() << "\n"
	   << "verifiedPartLnls="  << debugTraln->getPartitionLnls() << "\n"; 

#if HAVE_PLL != 0
      evaluateGeneric(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), traln.getAnyBranch().findNodePtr(traln), TRUE); 
#else 
      evaluateGeneric(&(traln.getTrHandle()),  traln.getAnyBranch().findNodePtr(traln), TRUE); 
#endif  

      tout << "after full traversal on orig=" << traln.getTrHandle().likelihood << std::endl; 

      // what to print? 
      assert(0);
      // tout << "current tree: " << traln << std::endl; 
      // tout << "help tree: " <<  *debugTraln << std::endl; 	        
    }  
  assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS); 
  debugPrint = 1;  
}
#ifdef DEBUG_LNL_VERIFY
void LikelihoodEvaluator::setDebugTraln(std::shared_ptr<TreeAln> _debugTraln)
{
  verifyLnl = true; 
  debugTraln = _debugTraln; 
}
#endif


void LikelihoodEvaluator::markDirty(const TreeAln &traln, nat nodeId)
{
  nat id = nodeId - traln.getNumberOfTaxa()  - 1;  
  for(nat i = 0; i < traln.getNumberOfPartitions() ;++i)
    arrayOrientation.setOrientation(i, id, INVALID); 
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
  // tout << "marking " << nodeId << " (" << id << ") dirty" <<  std::endl; 
  arrayOrientation.setOrientation(partitionId,id, INVALID); 
} 

void LikelihoodEvaluator::markPartitionDirty(const TreeAln &traln, nat partition )
{
  arrayOrientation.setPartitionInvalid(partition);
}

void LikelihoodEvaluator::accountForRejection(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat> &invalidNodes)
{
  arrayPolicy->accountForRejection(traln,partitions, invalidNodes, arrayOrientation); 
}

void LikelihoodEvaluator::freeMemory()
{ 
  arrayPolicy->freeMemory() ; 
#ifdef EVAL_DEBUG
  tout << arrayOrientation << std::endl; 
#endif
}  

void LikelihoodEvaluator::imprint(const TreeAln &traln) 
{ 
  arrayPolicy->imprint(traln , arrayOrientation);  
}
