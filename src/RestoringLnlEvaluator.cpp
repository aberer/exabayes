#include "RestoringLnlEvaluator.hpp"


RestoringLnlEvaluator::RestoringLnlEvaluator(std::shared_ptr<LnlRestorer> _restorer)
{
  this->restorer = _restorer; 
}


double RestoringLnlEvaluator::evaluate(TreeAln &traln, const Branch &evalBranch, bool fullTraversal )  
{
  int model  = ALL_MODELS; 
  nodeptr start = evalBranch.findNodePtr(traln) ; 

  restorer->traverseAndSwitchIfNecessary(traln, start, model, fullTraversal);
  restorer->traverseAndSwitchIfNecessary(traln, start->back, model, fullTraversal);
  
  exa_evaluateGeneric(traln,start,fullTraversal);   
#ifdef DEBUG_LNL_VERIFY
  if(verifyLnl)
    expensiveVerify(traln);
#endif
  return traln.getTr()->likelihood;  
}


// must be partial 
void RestoringLnlEvaluator::evalSubtree(TreeAln  &traln, const Branch &evalBranch)   
{ 
  bool masked = false; 
  nodeptr p = evalBranch.findNodePtr(traln); 
  int numberToExecute = 0; 
  int modelToEval = ALL_MODELS; 

  auto execModel = traln.getExecModel();
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      if(execModel[i])
	{
	  numberToExecute++;
	  modelToEval = i; 
	}
    }

  assert(numberToExecute == 1 || numberToExecute == traln.getNumberOfPartitions()); 
  if(numberToExecute > 1)
    modelToEval = ALL_MODELS; 

  if(p->x)
    {
      tree *tr = traln.getTr();
      assert(not isTip(p->number, tr->mxtips)); 
      p->x = 0; 
      p->next->x = 1; 
    }
  restorer->traverseAndSwitchIfNecessary(traln, p, modelToEval, false); 
  coreEvalSubTree(traln,p,masked); // NEEDED
}



double RestoringLnlEvaluator::evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions)    
{
  Branch root = findVirtualRoot(traln); 
  
  auto tr = traln.getTr(); 
  nat numPart = traln.getNumberOfPartitions(); 
  auto perPartitionLH = traln.getPartitionLnls();

  std::vector<bool> toExecute(numPart, false );   
  for(auto m : partitions)
    toExecute[m] = true; 
  traln.setExecModel(toExecute); 

  disorientTree(traln, root); 
  
  for(auto p : partitions)
    {
      restorer->traverseAndSwitchIfNecessary(traln,root.findNodePtr(traln), p, true ); 
      restorer->traverseAndSwitchIfNecessary(traln,root.getInverted().findNodePtr(traln), p, true); 
    }

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

#ifdef DEBUG_LNL_VERIFY
  exa_evaluateGeneric(traln, root.findNodePtr(traln),  FALSE );   
  if(fabs(prevLnl - traln.getTr()->likelihood) > 1e-6)
    {
      std::cout << "error while resetting lnl arrays. Likelihood should be " 
		<< prevLnl << " but was " << traln.getTr()->likelihood << std::endl; 
      assert(0); 
    }
#endif
}

void RestoringLnlEvaluator::imprint(const TreeAln &traln)
{ 
  prevLnl = traln.getTr()->likelihood; 
  restorer->resetRestorer(traln);   
}
