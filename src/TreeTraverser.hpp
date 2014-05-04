#ifndef _TREE_TRAVERSER_HPP
#define _TREE_TRAVERSER_HPP

#include <unordered_map>

#include "eval/LikelihoodEvaluator.hpp"
#include "model/TreeAln.hpp"
#include <functional>

class TreeTraverser
{
public: 
  TreeTraverser(bool doFirst, int depth, TreeAln &traln, LikelihoodEvaluator& eval, std::vector<AbstractParameter*>  params, const BranchPlain &rootOfTraversed, const BranchPlain &prunedSubtree); 

  std::unordered_map<BranchPlain, log_double> getResult() const {return _result; }
  void traverse() ; 

private : 			// METHODS
  void testInsert( const BranchPlain &insertBranch, BranchLengths floatingBranch,bool isFirst, int curDepth); 
  

private: 			// ATTRIBUTES
  bool _doFirst; 
  int _depth; 
  std::reference_wrapper<TreeAln> _traln; 
  std::reference_wrapper<LikelihoodEvaluator> _eval; 
  std::vector<AbstractParameter*> _params; 
  BranchPlain _rootOfTraversed; 
  BranchPlain _prunedSubtree; 
  
  std::unordered_map<BranchPlain, log_double> _result; 
  
  nat _cnt  ;			// for debug 
}; 


#endif

