#if 0 

#include "PrescoredDrawer.hpp"
#include "ParsimonyEvaluator.hpp"
#include "parameters/BranchLengthsParameter.hpp"




Branch2parsScores PrescoredDrawer::getPrescoredMap(TreeAln &traln ) 
{
  auto result = Branch2parsScores{};
  auto eval = ParsimonyEvaluator{}; 

  auto partitionPars = std::vector<nat>{}; 
  auto pLen = std::vector<nat>{}; 
  eval.evaluate(traln, traln.getAnyBranch().findNodePtr(traln), true, partitionPars, pLen);
  nat initScore  = std::accumulate(partitionPars.begin(),partitionPars.end(), 0); 
  tout << initScore << std::endl; 

  auto blParam = BranchLengthsParameter(0,0);
  auto param = &blParam; 
  auto params = std::vector<AbstractParameter*> {param};

  for(auto &branch : traln.extractBranches(param))
    {
      auto bvars = std::vector<Branch>{}; 
      if(not traln.isTipNode(branch.findNodePtr(traln)))
	bvars.push_back(branch) ;
      branch = branch.getInverted();
      if(not traln.isTipNode(branch.findNodePtr(traln)))
	bvars.push_back(branch) ;

      for(auto &b :  bvars)
	{
	  auto desc = traln.getDescendents(b); 

	  if(not desc.first.isTipBranch(traln))
	    {
	      
	    }
	}
      
    }
  
  return result; 
} 


Branch2Weight PrescoredDrawer::getWeights()
{
  auto result = Branch2Weight{}; 



  return result; 
}

#endif
