#ifndef TREE_INTEGRATOR
#define TREE_INTEGRATOR

#include <unordered_map>

#include "TreeAln.hpp"
#include "Chain.hpp"



typedef std::unordered_map<BranchLength,std::pair<double,double>> Branch2Stat; 

class TreeIntegrator
{
public: 
  TreeIntegrator(TreeAln& tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed);
  void prepareChain( const TreeAln &otherTree, bool copyEverything); 

  void integrateAllBranchesNew(const TreeAln &otherTree, std::string runid, nat sprDistance); 
  Branch2Stat integrateTree(double essThresh, LikelihoodEvaluator &eval ); 

private: 
  std::unique_ptr<Chain> integrationChain; 

}; 



#endif
