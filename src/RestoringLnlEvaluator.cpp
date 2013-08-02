#include "RestoringLnlEvaluator.hpp"

#include "GlobalVariables.hpp"


RestoringLnlEvaluator::RestoringLnlEvaluator(std::shared_ptr<ArrayRestorer> _restorer)  
{
  this->restorer = _restorer; 
}


double RestoringLnlEvaluator::evaluate(TreeAln &traln, const Branch &evalBranch, bool fullTraversal )  
{
  std::vector<nat> partitionsToEvaluate; 
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    partitionsToEvaluate.push_back(i); 

  nodeptr start = evalBranch.findNodePtr(traln) ; 

  if(fullTraversal)
    disorientTree(traln, evalBranch); 
  
  Branch root(start->number, start->back->number);
  restorer->toplevelSwitch(traln, root, partitionsToEvaluate, fullTraversal);

  exa_evaluateGeneric(traln,start,FALSE); // must be FALSE
#ifdef DEBUG_LNL_VERIFY
  if(verifyLnl)
    expensiveVerify(traln);
#endif
  return traln.getTr()->likelihood;  
}


// must be partial 
void RestoringLnlEvaluator::evalSubtree(TreeAln  &traln, const Branch &evalBranch)   
{ 
  assert(0); 

  // TODO does this do what you want? 
  
  bool masked = false; 
  nodeptr p = evalBranch.findNodePtr(traln); 
  // int numberToExecute = 0; 
  std::vector<nat> modelsToEval; 
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    modelsToEval.push_back(i); 

  if(p->x)
    {
      tree *tr = traln.getTr();
      assert(not isTip(p->number, tr->mxtips)); 
      p->x = 0; 
      p->next->x = 1; 
    }
  restorer->traverseAndSwitchIfNecessary(traln, p, modelsToEval, false); 
  coreEvalSubTree(traln,p,masked); 
}




double RestoringLnlEvaluator::evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions, bool fullTraversal)    
{
  // per-partition stuff not implemented yet (would be needed for per-partition bls)
  assert(fullTraversal); 

  Branch root = findVirtualRoot(traln); 
  
  auto tr = traln.getTr(); 
  nat numPart = traln.getNumberOfPartitions(); 
  auto perPartitionLH = traln.getPartitionLnls();

  std::vector<bool> toExecute(numPart, false );   
  for(auto m : partitions)
    toExecute[m] = true; 
  traln.setExecModel(toExecute); 

  disorientTree(traln, root); 
  
  restorer->toplevelSwitch(traln, root, partitions, fullTraversal);    

  exa_evaluateGeneric(traln, root.findNodePtr(traln),  FALSE ); 
  
  auto pLnl = traln.getPartitionLnls();
  for(auto m : partitions )
    perPartitionLH[m] = pLnl[m]; 

  traln.setPartitionLnls(perPartitionLH); 

  tr->likelihood = 0; 
  for_each(perPartitionLH.begin(), perPartitionLH.end(), [&](double d){ tr->likelihood += d; }); 
  traln.setExecModel(std::vector<bool>(numPart, true));

#ifdef DEBUG_LNL_VERIFY
  expensiveVerify(traln);   
#endif
  
  return tr->likelihood; 
}





void RestoringLnlEvaluator::resetToImprinted(TreeAln &traln)
{
  restorer->restoreArrays(traln); 
  Branch root = findVirtualRoot(traln); 

  auto p = root.findNodePtr(traln), 
    q = p->back; 

  assert( (traln.isTipNode(p) ||  p->x) 
	  && (traln.isTipNode(q) || q->x ) );

#if 0 
#ifdef DEBUG_LNL_VERIFY
  auto evalP = root.findNodePtr(traln); 

  disorientTree(traln, root) ;
  exa_evaluateGeneric(traln, evalP,  FALSE  );   

  traln.printArrayStart(tout);

  if(fabs(prevLnl - traln.getTr()->likelihood) > 1e-6)
    {                  
      std::cout << "error while resetting lnl arrays. Likelihood should be " 
  		<< prevLnl << " but was " << traln.getTr()->likelihood << std::endl; 

      // tout << traln << std::endl;
      for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
  	{
  	  auto p = traln.getPartition(i);
  	  tout << *p << std::endl; 
  	}
      assert(0); 
    }
#endif
#endif

  traln.getTr()->likelihood = prevLnl; 
  traln.setPartitionLnls(partitionLnls);   
}


void RestoringLnlEvaluator::imprint(const TreeAln &traln)
{ 
  prevLnl = traln.getTr()->likelihood; 
  partitionLnls = traln.getPartitionLnls(); 
  // tout << "imprinting lnl=" << prevLnl << std::endl; 
  restorer->resetRestorer(traln);   
}


void RestoringLnlEvaluator::resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions) 
{
  restorer->restoreSomePartitions(traln, partitions); 
}
