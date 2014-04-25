#include <array>

#include "system/BoundsChecker.hpp"
#include "LikelihoodSPR.hpp"
#include "TreeRandomizer.hpp"
#include "priors/AbstractPrior.hpp"
#include "GibbsProposal.hpp" 
#include "math/Density.hpp"


#define NUM_ITER 30   
#define NUM_BRANCH_OPT 3 
LikelihoodSPR::LikelihoodSPR(nat _minStep, nat _maxStep, double _likeWarp)
  : AbstractProposal(Category::TOPOLOGY, "likeSpr", 5., 0,0)
  ,  minStep(_minStep)
  , maxStep(_maxStep)
  , likeWarp(_likeWarp)
{
}


double LikelihoodSPR::scoreReattachment(TreeAln& traln, const BranchLength &reattachmentBranch, nodeptr prunedSubtree, LikelihoodEvaluator& eval, bool doRevert, bool isForward ) 
{ 
  double result = 0; 
  auto blParams = getBranchLengthsParameterView();
  assert(blParams.size() == 1 ); 
  auto blParam = blParams[0]; 

  auto q = reattachmentBranch.findNodePtr(traln), 
    r = reattachmentBranch.getInverted().findNodePtr(traln),
    p = prunedSubtree; 

  traln.clipNode(p->next,q); 
  traln.clipNode(p->next->next,r); 
  assert(p->next->back == q); 

  auto toOptimise = std::array<BranchLength,3>
    {
      {
	traln.getBranch(BranchPlain(p->number, p->back->number), blParam),
	traln.getBranch(BranchPlain(p->number, q->number), blParam), 
	traln.getBranch(BranchPlain(p->number, r->number),blParam)
      } 
    };  

  // EXAMINE: how many iter necessary? 
  auto mapPart = std::unordered_map<BranchLength, double>(); 
  for(nat i = 0;  i < NUM_BRANCH_OPT; ++ i)
    {
      // tout << MAX_SCI_PRECISION ; 
      // tout << "iter " << i << std::endl; 
      for(auto b : toOptimise)
	{
	  assert(blParams.size( )== 1 ); 
	  
	  auto result = GibbsProposal::optimiseBranch(traln,b, eval, NUM_ITER , blParams[0]);   
	  auto optBranch = result[0] ; 
	  auto nrd2 = result[2]; 
	  
	  // eval.evaluate(traln, b,true); 

	  // double optLen = optBranch.getInterpretedLength(traln, blParams[0]); 
	  // tout << curLen << "\t" << optLen << "\t" <<  PERC_PRECISION << ( optLen / curLen ) * 100   << MAX_SCI_PRECISION<< std::endl; 
	  traln.setBranch(optBranch, blParam); 
	  mapPart[optBranch] = nrd2; 
	}  
    }
  assert(mapPart.size() == toOptimise.size()); 

  map[reattachmentBranch] = mapPart; 

  // tout << std::endl; 
  eval.evaluate(traln, toOptimise[0].toPlain(), false); 
  result =  traln.getTrHandle().likelihood;
  
  // revert 
  if(doRevert)
    {      
      traln.detachNode(prunedSubtree); 
      traln.clipNode(q,r); 
      traln.setBranch(reattachmentBranch, blParam);
    }
  
  return result; 
}



BranchToLnlMap LikelihoodSPR::convertToProbMap(const BranchToLnlMap &map) const 
{
  auto result = map; 

  double best = std::numeric_limits<double>::lowest(); 
  for(auto &v : result)
    if( best < v.second ) 
      best = v.second; 

  // tout << MAX_SCI_PRECISION ; 

  // tout << std::endl << "after transform " << std::endl; 
  for(auto &v : result)
    {
      v.second -= best; 
      v.second = exp(v.second * likeWarp);
      // tout << v.second << std::endl; 
    }
  
  double sum = 0; 
  for(auto &v : result)
    sum += v.second; 
  
  if(std::isinf(sum))
    {
      tout << "problem with  summing up the following values: "  << std::endl; 
      for(auto &v : result)
	tout << v.second << std::endl; 
    }

  assert(not std::isinf(sum)); 
	 
  for(auto &v : result)
    v.second /= sum; 

  return result; 
}



std::pair<BranchLength,double> LikelihoodSPR::drawReattachment( const BranchToLnlMap &map, Randomness &rand) const 
{
  double r =  rand.drawRandDouble01(); 
  for( auto &v : map)
    {
      if(r < v.second)
	return v; 
      else 
	r -= v.second; 
    }

  tout << "we had " << map.size() << " options" << std::endl; 
  tout << "rest was " << r  << std::endl; 

  for(auto & v : map)
    {
      tout << v.first << "\t" << v.second << std::endl; 
    }
  
  assert(0); 
  return std::make_pair(BranchLength(0,0), -1);
}



BranchToLnlMap LikelihoodSPR::scoreOnBothSides(TreeAln &traln, const BranchLength &branch , nodeptr subtree, LikelihoodEvaluator &eval, bool isForward)  
{
  auto blParams = getBranchLengthsParameterView(); 
  assert(blParams.size( )== 1 ); 
  auto blParam = blParams[0]; 

  auto result = BranchToLnlMap(); 
  for(auto b : {branch, branch.getInverted()})
    {
      auto pTmp = b.findNodePtr(traln); 
      if(not traln.isTipNode(pTmp))
	{
	  auto desc = traln.getDescendents(b.toPlain());
	  scoreReattachmentInRadius(traln, traln.getBranch(desc.first, blParam).getInverted(), subtree, 1, eval, result, isForward); 
	  scoreReattachmentInRadius(traln, traln.getBranch(desc.second, blParam).getInverted(), subtree, 1, eval, result, isForward);       
	}
    }
  assert(result.size() > 1); 
  return result; 
}


BranchPlain LikelihoodSPR::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
{
  auto prunedSubtree = BranchPlain();
  nat depth = 0; 

  // determine pruned branch and neighbors 
  do
    {
      prunedSubtree = TreeRandomizer::drawBranchWithInnerNode(traln, rand);

      // // save brances adjacent to pruning point for restoration later. 
      // auto tmp  = traln.getDescendents(prunedSubtree.toPlain()); 
      // adjBranches.first = traln.getBranch(tmp.first, blParam).getInverted(); 
      // adjBranches.second = traln.getBranch(tmp.second, blParam).getInverted(); 
      
      // p = prunedSubtree.findNodePtr(traln);
      // q = adjBranches.first.findNodePtr(traln); 
      // r = adjBranches.second.findNodePtr(traln); 

      depth = traln.getDepth(prunedSubtree.getInverted().toPlain()); 
    } while( int(depth)-2 < int(minStep)  );
  // TODO (really) think the stop condition through again 

  assert(depth > 2 ); 

  return prunedSubtree; 
} 


void LikelihoodSPR::determineMove(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
  // tout << "tree before: " << TreePrinter(false, true, false).printTree(traln) << std::endl; 
  // tout << MAX_SCI_PRECISION << "initLnl=" << traln.getTr()->likelihood << std::endl;  

  auto blParams = getBranchLengthsParameterView();   
  assert(blParams.size( )== 1 ); 
  auto blParam = blParams[0]; 

  auto prunedSubtreeTmp = determinePrimeBranch(traln, rand); 
  auto prunedSubtree = traln.getBranch(prunedSubtreeTmp, blParam); 

  // save brances adjacent to pruning point for restoration later. 
  auto tmp  = traln.getDescendents(prunedSubtree.toPlain()); 
  auto adjBranches = std::make_pair(BranchLength(0,0), BranchLength(0,0)); 
  adjBranches.first = traln.getBranch(tmp.first, blParam).getInverted(); 
  adjBranches.second = traln.getBranch(tmp.second, blParam).getInverted(); 
      
  auto p = prunedSubtree.findNodePtr(traln);
  auto q = adjBranches.first.findNodePtr(traln); 
  auto r = adjBranches.second.findNodePtr(traln); 

  // auto q = nodeptr(nullptr); 
  // auto r = nodeptr(nullptr); 
  // auto p = nodeptr(nullptr);

  // tout << "pruning " << prunedSubtree <<  " from " << adjBranches.first << "\t" << adjBranches.second << std::endl; 
  // tout << "depth below  " << prunedSubtree << " is " << depth << std::endl; 
  // tout << "path: " << traln.getLongestPathBelowBranch(prunedSubtree.getInverted()) << std::endl; 

  
  // initial prune 
  traln.clipNode(q,r);
  traln.detachNode(p); 
  auto lengthsA = adjBranches.first.getLength(); 
  auto lengthsB = adjBranches.second.getLength(); 
  // auto lengths = std::vector<double>();  
  // for(nat i = 0; i < lengthsA.size(); ++i)    
    // lengths.push_back(lengthsA.at(i) * lengthsB.at(i)); 
  auto length = lengthsA * lengthsB; 
  auto branchAfterPrune = BranchLength(q->number, r->number );
  branchAfterPrune.setLength(length);
  traln.setBranch(branchAfterPrune, blParam);
  
  // somehow treat the branch after pruning  
  {

    assert(blParams.size() == 1); 
    auto blParam =  blParams[0]; 
    auto result = GibbsProposal::optimiseBranch( traln, branchAfterPrune, eval, 30,  blParam); 
    auto optBranch = result[0]; 
    traln.setBranch(optBranch, blParam) ;
  }

  // get relevant insertion positions on both sides of the pruning
  // branch
  auto forwardMap = scoreOnBothSides(traln, branchAfterPrune, p, eval, true );

  // remove initial branch (we want to draw a new one)
  assert(forwardMap.find(branchAfterPrune) == forwardMap.end( )); 
  auto forwardProbs = convertToProbMap(forwardMap); 

  // for(auto &elem : forwardMap)
  //   tout << elem.first <<  "\t" << elem.second << "\t" << forwardProbs[elem.first] << std::endl;         
  auto drawnElem = drawReattachment(forwardProbs, rand);
  auto drawnPos = traln.getBranch(drawnElem.first.toPlain(),blParams); 

  // tout << drawnPos  << std::endl;   

  // optimise for the backward move 
  {
    // auto optBranch = GibbsProposal::optimiseBranch(traln,drawnPos, eval, NUM_ITER, blParams); 
    // double nrd1 = 0, nrd2 = 0; 

    // HACK 
    auto bl = traln.getBranch(drawnPos.toPlain(),blParam); 

    auto result = GibbsProposal::optimiseBranch( traln, bl, eval, NUM_ITER,  blParam); 
    auto optBranch = result[0]; 
    // , nrd1, nrd2
    
    traln.setBranch(optBranch, blParam) ;
  }

  // insert into the  drawn position 
  // important TODO: make cheaper: no evaluation? save optimized branches  
  auto backwardMap = scoreOnBothSides(traln, traln.getBranch(drawnPos.toPlain(),blParam), p, eval, false);
  assert(backwardMap.find(branchAfterPrune) != backwardMap.end()); 
  auto backwardProb = convertToProbMap(backwardMap); 
  auto backProb = backwardProb.at(branchAfterPrune); 

  // tout << "backP=" << log(backProb) << "\tforP="<< log(drawnElem.second) << std::endl; 

  hastings += log(backProb) - log(drawnElem.second) ; 

  // tout << "hastings="<< hastings <<  std::endl; 

  traln.setBranch(drawnPos, blParams);
  
  // restore 
  traln.clipNode(p->next, r); 
  traln.clipNode(p->next->next, q); 

  assert(blParams.size() == 1 ); 

  traln.setBranch(adjBranches.first, blParam); 
  traln.setBranch(adjBranches.second, blParam); 
  traln.setBranch(prunedSubtree, blParam); 

  // tout  << "tree after: " << traln << std::endl; ; 
  // eval.evaluate(traln, prunedSubtree, treu); 
  // tout << "extracting from " << prunedSubtree << "\t" << drawnPos << std::endl; 
  auto bls = std::make_tuple(prunedSubtree.toPlain(), drawnPos.toPlain() ); 
  move.extractMoveInfo(traln, bls , blParams);
  // tout << "determined move: " << move << std::endl; 
  // tout << "lnl after: " << eval.evaluate(traln, prunedSubtree, true) << std::endl;

  // TODO prior 
  // tout << std::endl; 
}



void LikelihoodSPR::proposeBranches(TreeAln& traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval )
{
#if 0 
  auto blParams = getBranchLengthsParameterView(); 
  auto& relMap = map.at(move.getInsertionBranchBefore().toBlDummy());

  // eval.evaluate(traln, move.getEvalBranch(traln), false);

  assert(blParams.size( )== 1 ); 
  auto param = blParams[0]; 

  tout << MAX_SCI_PRECISION; 

  // handle the outer pruning branch first  
  {
    auto afterPrune = traln.getBranch(move.getPruningBranchAfterPrune(), param); 

    // double beforeInternal = afterPrune.getLength(); 
    double prevLen = afterPrune.getInterpretedLength(traln, param);
    double hastingsPart = 0; 
    auto optBranch = GibbsProposal::drawFromEsitmatedPosterior(afterPrune, eval, traln, rand, NUM_ITER, hastingsPart, param);
    hastings += hastingsPart;  


    traln.setBranch(optBranch, param); 
    // double afterInternal = optBranch.getLength(); 
    double curLen = optBranch.getInterpretedLength(traln,param); 
    
    prior.addToRatio( 
		     param->getPrior()->getLogProb( { curLen } )
		     - param->getPrior()->getLogProb( { prevLen } )
		      ); 
  }

  // first one is the branch existing after move 
  auto branches = std::array<BranchPlain,3> 
    { {
      move.getSubtreeBranchAfter(traln), // std::make_pair(, move.getSubtreeBranch(traln)),
      move.getInsertionBranchAfterInner(), // std::make_pair(move.getPruningBranchBeforeInner(), ),
      move.getInsertionBranchAfterOuter() // std::make_pair(move.getInsertionBranchBefore(), )      
      }} ; 
  
  // propose for these branches 
  for(auto branch : branches)
    {
      auto bl = traln.getBranch(branch, param ); 
      auto lengthBefore = bl.getInterpretedLength(traln,param); 
      auto iter = relMap.find(bl); 
      double nropt = iter->first.getInterpretedLength(traln, param);
      // double afterInternal = iter->first.getLength(param); 
      double nrd2 = iter->second; 

      auto proposalResult = GibbsProposal::propose(nropt, nrd2, rand); 
      double proposal = proposalResult[0]; 
      double alpha = proposalResult[1]; 
      double beta = proposalResult[2]; 

      auto newBranch = branch.toBlDummy();  
      newBranch.setConvertedInternalLength(traln, param, proposal);
      if(not BoundsChecker::checkBranch(newBranch))
	BoundsChecker::correctBranch(newBranch); 
      traln.setBranch(newBranch, param); 


      auto prevLen = bl.getInterpretedLength(traln, param); 
      auto curLen = newBranch.getInterpretedLength(traln, param); 

      auto pr = param->getPrior(); 
      prior.addToRatio( pr-> getLogProb( { curLen} )  -  pr->getLogProb( { prevLen} ) ) ; 
      
      // prior.updateBranchLengthPrior(traln, beforeInternal,newBranch.getLength(), param);
      
      double lnBackP = Density::lnGamma(lengthBefore, alpha, beta); 
      double lnForP = Density::lnGamma(proposal, alpha, beta); 

      double hastingsHere = lnBackP - lnForP; 
      // tout << "lengthBefore=" << lengthBefore << "\tproposal="  << proposal << "\thastings=" << hastingsHere << std::endl; 
      hastings += hastingsHere; 
    }
  
  // move.disorientAtNode(traln, move.getEvalBranch(traln).findNodePtr(traln));
  // eval.evaluate(traln, move.getEvalBranch(traln), false);
#endif
}



void LikelihoodSPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto multiplyBranchesUsingPosterior = 
#ifdef PROPOSE_BRANCHES_FOR_SPR 
    true; 
#else 
  false ; 
#endif


  auto blParams = getBranchLengthsParameterView();
  map.clear(); 

  determineMove(traln, prior, hastings, rand, eval);
  
  if(multiplyBranchesUsingPosterior)
    {  
      // NOTICE: extended! 
      // move.proposeBranches(traln, blParams, eval, hastings, rand, false); 
    }

  // auto blParam = blParams[0]; 
  // assert(blParams.size() == 1); 

#if 0 
  savedSubtreeBranch = traln.getBranch(move.getSubtreeBranchBefore(traln),blParam); 
#endif
  
  // tout << "nni-dist=" << move.getNniDistance() << std::endl; 
  move.applyToTree(traln, blParams);
  

  assert(0); 
  // why did we disorient here? 
  // for(auto &elem : move.get)
  // move.disorientAtNode(traln,move.getSubtreeBranchAfter(traln).findNodePtr(traln));

#if  0 
  proposeBranches(traln, prior, hastings, rand, eval);
#else 
  auto params = getBranchLengthsParameterView();
  // assert(params.size() == 1 ); 
  // auto param = params[0]; 
    // auto proposedBranches = move.proposeBranches(traln, blParams, eval, hastings, rand, true   );
  // for(auto &b : proposedBranches)
  //   {
  //     for(auto &param : params )
  // 	{
  // 	  auto prevLen = traln.getBranch(b.toPlain(), param).getInterpretedLength(traln , param); 
  // 	  auto dummy = BranchLength(); 
  // 	  dummy.setLength(b.getLength(param));
  // 	  auto curLen = dummy.getInterpretedLength(traln, param); 
	  
  // 	  auto pr = param->getPrior(); 
  // 	  auto ratio =  pr->getLogProb( { curLen } ) - pr->getLogProb( { prevLen } ) ; 
  // 	  prior.addToRatio(ratio );
	  
  // 	  // prior.updateBranchLengthPrior(traln, traln.getBranch(b.toPlain(), param).getLength(), b.getLength(param), param); 
  // 	}
  //     traln.setBranch(b,params);
  //   }
#endif
}



void  LikelihoodSPR::scoreReattachmentInRadius(TreeAln &traln, BranchLength attachment, nodeptr prunedSubtree, nat depth, LikelihoodEvaluator& eval, BranchToLnlMap &map, bool isForward) 
{
  auto blParams = getBranchLengthsParameterView(); 

  // for(nat i = 0; i < depth; ++i)
  //   tout << " "; 

  if( depth == 0  || (minStep <= depth && depth <= maxStep )  )
    {
      auto score = scoreReattachment(traln, attachment,  prunedSubtree, eval, true, isForward ); 
      assert(map.find(attachment) == map.end()); 
      map[attachment] =  score;       
      // tout << "inserting " <<  attachment << "\t" << score << std::endl; 
    }
  else 
    {
      // tout << "NOT inserting " << attachment << std::endl; 
    }
  
  if( not traln.isTipNode(attachment.findNodePtr(traln)) )
    {
      assert( blParams.size() == 1 ) ;
      auto blParam = blParams[0];

      auto desc =  traln.getDescendents(attachment.toPlain() );
      auto blA = traln.getBranch(desc.first.getInverted(), blParam); 
      auto blB = traln.getBranch(desc.second.getInverted(), blParam); 

      if(depth <=  maxStep)
	{
	  scoreReattachmentInRadius(traln, blA, prunedSubtree, depth+1, eval, map, isForward);
	  scoreReattachmentInRadius(traln, blB, prunedSubtree, depth+1, eval, map, isForward);
	}
    }
}


void LikelihoodSPR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  for(auto &elem : move.getDirtyNodes(traln, false))
    evaluator.markDirty( traln, elem);

  auto evalBranch = traln.getAnyBranch(); 

  // TODO inefficient 
  evaluator.evaluate(traln, evalBranch, true);
} 


void LikelihoodSPR::resetState(TreeAln &traln)  
{
  auto blParams = getBranchLengthsParameterView(); 
  auto blParam = blParams[0]; 
  assert(blParams.size() == 1); 
  move.revertTree(traln, blParams);
  traln.setBranch(savedSubtreeBranch, blParam);
} 


void LikelihoodSPR::writeToCheckpointCore(std::ostream &out)const 
{
  // noop 
}
 

void LikelihoodSPR::readFromCheckpointCore(std::istream &in) 
{
  // noop 
} 


AbstractProposal* LikelihoodSPR::clone() const 
{
  return new LikelihoodSPR(*this);
}



std::vector<nat> LikelihoodSPR::getInvalidatedNodes(const TreeAln& traln) const
{
  return move.getDirtyNodes(traln, false);
} 
