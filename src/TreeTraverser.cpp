#include "TreeTraverser.hpp"

#define PRINT_INFO 



TreeTraverser::TreeTraverser(bool doFirst, int depth, TreeAln &traln, LikelihoodEvaluator& eval, std::vector<AbstractParameter*>  params, const BranchPlain &rootOfTraversed, const BranchPlain &prunedSubtree)
  : _doFirst{doFirst}
  , _depth(depth)
  , _traln(traln)
  , _eval(eval)
  , _params{params}
  , _rootOfTraversed{rootOfTraversed}
  , _prunedSubtree{prunedSubtree}
  , _cnt{1}
{
#ifdef PRINT_INFO
  tout << "\n\ncreated tree traverser for " << SHOW(prunedSubtree) << SHOW(rootOfTraversed) << std::endl; 
#endif
}


void TreeTraverser::traverse()
{
  auto floatingBranch = BranchLengths(0,0); 
  auto branchAfterPrune = BranchPlain(0,0);

  std::tie(floatingBranch, branchAfterPrune) = _traln.get().pruneSubtree(_prunedSubtree, _rootOfTraversed, _params); 

// #ifdef PRINT_INFO
//   tout << "\n\ntestInsert starts with "<< branchAfterPrune << std::endl; 
// #endif
  
  testInsert(branchAfterPrune, floatingBranch, true, 0 ); 

  _traln.get().insertSubtree( _prunedSubtree, branchAfterPrune, floatingBranch, _params );
}


/** 
    insert, do stuff and remove again 
 */ 
void TreeTraverser::testInsert( const BranchPlain &insertBranch, BranchLengths floatingBranch, bool isFirst, int curDepth)
{
  if(_depth < curDepth )
    return; 

  if(not isFirst || _doFirst )
    {
      // modify the floating branch, s.t it tells where to insert.  We
      // assume that primNode of insertBranch is the direction we are
      // coming from (while descending), so this is the point, where
      // we have to insert the branch length

      floatingBranch.setPrimNode(_prunedSubtree.getPrimNode()); 
      floatingBranch.setSecNode( insertBranch.getPrimNode() ); 
      
      _traln.get().insertSubtree( _prunedSubtree,  insertBranch, floatingBranch, _params); 

      auto dirtyNodes = std::vector<nat>{ _prunedSubtree.getPrimNode() ,  insertBranch.getPrimNode()}; 
      tout << "dirty nodes: " << dirtyNodes << std::endl; 

      for(auto v : dirtyNodes)
	_eval.get().markDirty(_traln.get(), v); 

      _eval.get().evaluate(_traln.get(), _prunedSubtree, false);
      
      tout  << _cnt++ << "================================================================"  << std::endl; 

      // should not have this yet 
      auto lnl = _traln.get().getLikelihood(); 
      assert(_result.find(insertBranch) == _result.end()); 
      _result.insert(std::make_pair(insertBranch, lnl));

      auto floatAfterPrune = BranchLengths(0,0); 
      auto branchAfterPrune = BranchPlain(0,0); 

      std::tie( floatAfterPrune, branchAfterPrune) = _traln.get().pruneSubtree(_prunedSubtree, floatingBranch.toPlain(), _params); 
      
      // assure correct pruning 
      auto lenA = floatAfterPrune.getLengths(); 
      auto lenB = floatingBranch.getLengths(); 
      for(auto i = 0u; i < lenA.size() ; ++i)
	assert(lenA.at(i) == lenB.at(i)); 
    }

  if(not  _traln.get().isTipNode(insertBranch.getSecNode())  ) // insertBranch.isTipBranch(_traln.get())
    {
      auto desc = _traln.get().getDescendents(insertBranch.getInverted()); 

      testInsert(desc.first, floatingBranch, false, curDepth+1 ); 
      testInsert(desc.second, floatingBranch, false, curDepth+1 ); 
    }
}
