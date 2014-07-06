#ifndef TREE_INTEGRATOR
#define TREE_INTEGRATOR

#include <unordered_map>

#include "TreeAln.hpp"
#include "Chain.hpp"

typedef std::tuple<BranchPlain,double,double> BranchStat   ; 

class TreeIntegrator
{
public: 
  TreeIntegrator(TreeAln& traln, std::shared_ptr<TreeAln> debugTree, randCtr_t seed, ParallelSetup* plPtr); 
  void prepareChain( const TreeAln &otherTree, bool copyEverything); 
  std::vector<BranchStat> integrateTree(double essThresh, LikelihoodEvaluator &eval, std::vector<BranchPlain> branches); 
  void integrateAllMoves(const TreeAln &otherTree, std::string runid, nat sprDistance); 


private: 
  std::unique_ptr<Chain> integrationChain; 

}; 

#endif
