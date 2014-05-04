#include "LikelihoodSPR.hpp" 

#include "TreeTraverser.hpp"

#define PRINT_INFO

#define MY_EPS  1e-100

LikelihoodSPR::LikelihoodSPR(nat maxStep, double likeWarp)
  : AbstractProposal{Category::TOPOLOGY, "likeSPR", 5.,0.,0.,false}
  ,  _maxStep{maxStep}
  , _likeWarp{likeWarp}
{
} 


auto LikelihoodSPR::transformLikelihoods(std::unordered_map< BranchPlain, log_double > branch2lnl ) const 
-> std::unordered_map<BranchPlain,double> 
{
  auto result = std::unordered_map<BranchPlain,double>{} ; 

  auto lnlPairs = std::vector< std::pair<BranchPlain,log_double> >(begin(branch2lnl), end(branch2lnl));
  auto tmp = std::vector<log_double>(lnlPairs.size());
  std::transform( begin(lnlPairs), end(lnlPairs)  , begin(tmp), []( const std::pair<BranchPlain,log_double> &elem ) { return elem.second ; } ) ; 
  auto maximum = *max_element(begin(tmp), end(tmp));
  
  for(auto &v : lnlPairs) 
    v.second /= maximum; 

  for(auto & v : lnlPairs)
    v.second = exponentiate(v.second, _likeWarp);

  auto sum = 0.; 
  for(auto v : lnlPairs)
    {
      double normVal =  v.second.toAbs() + MY_EPS; 
      result.insert(std::make_pair(v.first, normVal)); 
      sum += normVal; 
    }

  for(auto &v : result)
    v.second /= sum; 

#ifdef PRINT_INFO
  tout << "================================================================" << std::endl; 
  tout << "likelihoods: " << std::endl; 
  for(auto &pair : result)
    tout << pair.first << "\t" <<  MAX_SCI_PRECISION << pair.second << "\t" << SOME_FIXED_PRECISION << pair.second << std::endl; 
#endif

  return result; 
}



auto  LikelihoodSPR::determineSprMove(TreeAln &traln, LikelihoodEvaluator &eval, Randomness &rand, const BranchPlain& prunedTree, std::vector<AbstractParameter*> &params)  
  -> std::tuple<SprMove, double>
{
  auto result = SprMove();
  auto branch2lnl = computeLikelihoodsOfInsertions(traln,eval, prunedTree, params);
  auto branch2prob = transformLikelihoods(branch2lnl);

  auto r = rand.drawRandDouble01(); 

  auto choice = []( const std::unordered_map<BranchPlain, double> &map, double r ) -> BranchPlain
    {
      for(auto &v : map)
	{
	  r -= v.second; 
	  if(r <= 0 )
	    return v.first; 
	}

      assert(0); 
      return BranchPlain(0,0); 
    }; 
  
  auto chosenBranch = choice(branch2prob, r); 

  auto prob = branch2prob[chosenBranch]; 

#ifdef PRINT_INFO  
  tout << SHOW(chosenBranch)  << std::endl;
#endif
  
  result.extractMoveInfo(traln, std::make_tuple(prunedTree, chosenBranch), params ); 
  
  return std::make_pair(result, prob); 
} 


auto  LikelihoodSPR::computeLikelihoodsOfInsertions(TreeAln &traln, LikelihoodEvaluator &eval, const BranchPlain& prunedTree, std::vector<AbstractParameter*> &params)  
  -> std::unordered_map<BranchPlain, log_double>
{
  eval.evaluateSubtrees(traln,prunedTree,false);		      // TODO efficiency?

  auto desc = traln.getDescendents(prunedTree); 

  auto resultA = std::unordered_map<BranchPlain,log_double>{}; 
  auto resultB = std::unordered_map<BranchPlain,log_double>{}; 

  assert( not desc.first.isTipBranch(traln) ||  not desc.second.isTipBranch(traln)); 
  if(not desc.first.isTipBranch(traln))
    {
      auto tt = TreeTraverser( false, _maxStep, traln, eval, params, desc.first, prunedTree ); // TODO param: branches already known! 
      tt.traverse();
      resultA = tt.getResult(); 
    }

  if(not desc.second.isTipBranch(traln))
    {
      auto tt = TreeTraverser(false,  _maxStep, traln, eval, params, desc.second, prunedTree);
      tt.traverse();
      resultB = tt.getResult(); 
    }


  auto expSize = resultA.size() + resultB.size(); 
  // join the maps 
  resultA.insert(begin(resultB), end(resultB)); 
  assert(expSize == resultA.size()); 

  return resultA;  
}




void LikelihoodSPR::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto params = getBranchLengthsParameterView(); 
  auto prunedTree = determinePrimeBranch(traln, rand); 

#ifdef PRINT_INFO  
  tout << SHOW(prunedTree) << std::endl; 
#endif
  
  auto forwProb = 0.; 
  std::tie( _move, forwProb) = determineSprMove(traln, eval, rand, prunedTree, params);
  
#ifdef PRINT_INFO
  tout <<  SHOW(_move) << std::endl; 
#endif
  tout << "================================================================" << std::endl; 
  tout << "================================================================" << std::endl; 

  // maybe do branch length proposal  

  _move.applyToTree(traln, params);
  auto nodes = _move.getDirtyNodes(traln, false); 
  for(auto v : nodes)
    eval.markDirty(traln, v);
  
  auto lnlMapBack = computeLikelihoodsOfInsertions(traln, eval, prunedTree, params); 
  auto probMapBack = transformLikelihoods(lnlMapBack); 

  auto invMove = _move.getInverseMove(traln, params); 
  auto backBranch =  invMove.getInsertBranch(); 
  tout << SHOW(backBranch) << std::endl; 

  // tout << "back insertions: " << std::endl; 
  // for(auto& v : probMapBack)
  //   tout << v.first << "\t" << MAX_SCI_PRECISION << v.second << "\t"<< SOME_FIXED_PRECISION << v.second  << std::endl; 
  
  assert(probMapBack.find(backBranch) != end(probMapBack)); 

  auto iter = probMapBack.find(backBranch); 
  auto backProb = iter->second; 

  // tout << MAX_SCI_PRECISION << SHOW(backProb) << SHOW(forwProb)<< std::endl ; 
  hastings *= log_double::fromAbs(backProb / forwProb); 
} 


void LikelihoodSPR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  auto toEval = _move.getEvalBranch(traln);

  for(auto &elem : _move.getDirtyNodes(traln, false))
    evaluator.markDirty(traln, elem); 

  evaluator.evaluate(traln,toEval, false); 
} 


void LikelihoodSPR::resetState(TreeAln &traln)   
{
  auto params = getBranchLengthsParameterView();
  _move.revertTree(traln, getSecondaryParameterView() ); 
}


BranchPlain LikelihoodSPR::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
{
  auto prunedTree = BranchPlain();
  nodeptr p, pn, pnn;  
  do 
    {
      prunedTree  = TreeRandomizer::drawBranchWithInnerNode(traln,rand); 
      p = prunedTree.findNodePtr(traln);
      pn = p->next->back; 
      pnn = p->next->next->back;         

    } while( (traln.isTipNode(pn) &&  traln.isTipNode(pnn))     ); 

  return prunedTree; 
} 


auto LikelihoodSPR::getInvalidatedNodes(const TreeAln &traln) const 
-> std::vector<nat>
{
  return _move.getDirtyNodes(traln, false); 
}  



void LikelihoodSPR::printParams(std::ostream &out)  const 
{
  out << " ; radius=" << _maxStep << "; warp="<< _likeWarp ; 
}
