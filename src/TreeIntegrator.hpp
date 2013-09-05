#ifndef TREE_INTEGRATOR
#define TREE_INTEGRATOR

#include <unordered_map>

#include "TreeAln.hpp"
#include "Chain.hpp"



typedef std::unordered_map<BranchLength,std::pair<double,double>> Branch2Stat; 

class TreeIntegrator
{
public: 
  TreeIntegrator(std::shared_ptr<TreeAln> tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed);
  void prepareChain( const TreeAln &otherTree, bool copyEverything); 

  void integrateAllBranches(const TreeAln &traln, std::string runid);   
  Branch2Stat integrateTree(nat numgen, nat thinning); 
   
private: 
  std::unique_ptr<Chain> integrationChain; 

}; 



#endif
