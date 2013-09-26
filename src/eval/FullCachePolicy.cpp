#include "FullCachePolicy.hpp"


FullCachePolicy::FullCachePolicy(const TreeAln& traln, bool cacheTipTip, bool cacheTipInner  )
  : restorer(traln, cacheTipTip, cacheTipInner)
{
}


void FullCachePolicy::imprintPolicy(const TreeAln &traln, ArrayOrientation &arrayOrient)  
{
  restorer.resetRestorer(traln, arrayOrient);   
}


void FullCachePolicy::freeMemory()   
{
  restorer.clearMemory();
}

void FullCachePolicy::accountForRejectionPolicy(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat>& invalidNodes, ArrayOrientation &arrayOrient)
{
  restorer.restoreSomePartitions(traln, partitions, arrayOrient); 
}


std::unique_ptr<ArrayPolicy> FullCachePolicy::clone() const  
{
  return std::unique_ptr<ArrayPolicy>(new FullCachePolicy(*this) ); 
}

void FullCachePolicy::prepareForEvaluation(TreeAln &traln, BranchPlain virtualRoot, nat model, ArrayOrientation& arrayOrientation )
{
  restorer.traverseAndCache(traln, virtualRoot.findNodePtr(traln), model,  arrayOrientation );
} 
