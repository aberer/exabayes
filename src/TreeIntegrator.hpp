#ifndef TREE_INTEGRATOR
#define TREE_INTEGRATOR

#include <unordered_map>

#include "model/TreeAln.hpp"
#include "mcmc/Chain.hpp"



typedef std::unordered_map<BranchLength,std::pair<double,double>> Branch2Stat; 

class TreeIntegrator
{
public: 
  TreeIntegrator(TreeAln& traln, std::shared_ptr<TreeAln> debugTree, randCtr_t seed, ParallelSetup* plPtr); 
  // TreeIntegrator(TreeAln& tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed, std::shared_ptr<ParallelSetup> plPtr);
  void prepareChain( const TreeAln &otherTree, bool copyEverything); 

  void integrateAllBranchesNew(const TreeAln &otherTree, std::string runid, nat sprDistance); 
  Branch2Stat integrateTree(double essThresh, LikelihoodEvaluator &eval ); 

private: 
  std::unique_ptr<Chain> integrationChain; 

}; 



#endif
