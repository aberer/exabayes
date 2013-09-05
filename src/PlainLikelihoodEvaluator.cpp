#include "PlainLikelihoodEvaluator.hpp"

#include "Branch.hpp"

#include <cassert>


double PlainLikelihoodEvaluator::evaluate(TreeAln &traln, const BranchPlain &evalBranch,  bool fullTraversal )  
{
  auto evalP = evalBranch.findNodePtr(traln);

  exa_evaluateGeneric(traln,evalP,TRUE );   
#ifdef DEBUG_LNL_VERIFY
  if(verifyLnl)
    expensiveVerify(traln);
#endif

  return traln.getTr()->likelihood; 
}


// double PlainLikelihoodEvaluator::evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions)
// {
//   Branch root = findVirtualRoot(traln); 
  
//   auto tr = traln.getTr(); 
//   nat numPart = traln.getNumberOfPartitions(); 
//   auto perPartitionLH = traln.getPartitionLnls();

//   std::vector<bool> toExecute(numPart, false );   
//   for(auto m : partitions)
//     toExecute[m] = true; 
//   traln.setExecModel(toExecute); 

//   disorientTree(traln, root); 

//   exa_evaluateGeneric(traln, root.findNodePtr(traln),  FALSE ); 
  
//   auto pLnl = traln.getPartitionLnls();
//   for(auto m : partitions )
//     perPartitionLH[m] = pLnl[m]; 

//   traln.setPartitionLnls(perPartitionLH); 

//   tr->likelihood = 0; 
//   for_each(perPartitionLH.begin(), perPartitionLH.end(), [&](double d){ tr->likelihood += d; }); 
//   traln.setExecModel(std::vector<bool>(numPart, true));

// #ifdef DEBUG_LNL_VERIFY
//   expensiveVerify(traln);   
// #endif
  
//   return tr->likelihood; 
// }
 
void PlainLikelihoodEvaluator::evalSubtree( TreeAln &traln, const BranchPlain &evalBranch)
{ 
  // meh 
  assert(0); 

#if 0 
  bool masked = false; 
  nodeptr p = evalBranch.findNodePtr(traln); 
  int numberToExecute = 0; 

  auto execModel = traln.getExecModel();
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      if(execModel[i])
	numberToExecute++;
    }

  assert(numberToExecute == 1 || numberToExecute == traln.getNumberOfPartitions()); 

  if(p->x)
    {
      tree *tr = traln.getTr();
      assert(not isTip(p->number, tr->mxtips)); 
      p->x = 0; 
      p->next->x = 1; 
    }

  coreEvalSubTree(traln,p,masked); // NEEDED
#endif
} 


bool PlainLikelihoodEvaluator::traverseAndTouch(const TreeAln &traln, const BranchPlain &b)
{
  if(b.isTipBranch(traln))
    return true; 
  
  auto p = b.findNodePtr(traln); 
  bool result =  ( p->x == 1  ) ; 
  
  auto descending = traln.getDescendents(b); 
  result &= traverseAndTouch(traln,descending.first);
  result &= traverseAndTouch(traln,descending.second);

  if(not result )
    disorientNode(p);

  return result; 
} 



void PlainLikelihoodEvaluator::resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions) 
{
  assert(0); 
  // TODO 
}
