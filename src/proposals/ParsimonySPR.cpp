#include <unordered_map>	
#include <functional>

#include "priors/AbstractPrior.hpp"

#include "ParsimonySPR.hpp"
#include "model/Branch.hpp"

#include "comm/ParallelSetup.hpp"

#include "TreePrinter.hpp"

#include "BranchLengthOptimizer.hpp"
#include "GammaProposer.hpp"

static double state2factor(nat states)
{
  static double typicalBL = 0.05; 
  return log((1.0/states) - exp(-(states/(states-1) * typicalBL)) / states);
}


std::array<double,2> ParsimonySPR::factors = 
  {
    {
      state2factor(4),		// DNA
      state2factor(20)		// AA
    }
  }; 


// #define PRINT_DEBUG_PARS

ParsimonySPR::ParsimonySPR( double parsWarp, double blMulti, int depth)
  : AbstractProposal(Category::TOPOLOGY, "parsSPR", 5. ,  0,0, false)
  , _parsWarp(parsWarp)    
  , _blMulti(blMulti)    
  , _depth(depth)
{
  // very relevant to efficient, setting this as a varaible, s.t. we do not get lo t
  _parallelReduceAtEnd = true; 

  auto res = getEnvironmentVariable("PARS_BL"); 
  _proposePrimeBranch = res.compare("1") == 0 ; 
}


void ParsimonySPR::testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, branch22states2score &result, int curDepth,  const branch22states2score& alreadyComputed)
{
  if(curDepth == 0 )
    return; 
  --curDepth; 

  nodeptr insertBack =  insertPos->back;   
  traln.clipNode(insertPos, prunedTree->next);
  traln.clipNode( insertBack, prunedTree->next->next); 
  
  auto b = BranchPlain(insertPos->number, insertBack->number); 

  if(alreadyComputed.size() > 0 && alreadyComputed.find(b) != alreadyComputed.end()) 
    {
      result[b] = alreadyComputed.at(b); 
    }
  else 
    {
      ParsimonyEvaluator::disorientNode(prunedTree); 
      auto states2parsimony =   _pEval.evaluate(traln, prunedTree, false  ); 
      // assert(result.find(b) == result.end()) ;   
      result[b] = states2parsimony; 
    }

  traln.clipNode(insertPos, insertBack); 
  prunedTree->next->back = prunedTree->next->next->back = NULL; 

  // recursively descend 
  if(not traln.isTipNode(insertPos))
    {
      testInsertParsimony(traln, insertPos->next->back, prunedTree, result, curDepth, alreadyComputed); 
      testInsertParsimony(traln, insertPos->next->next->back, prunedTree, result, curDepth, alreadyComputed); 
    }
}


weightMap ParsimonySPR::getWeights(const TreeAln& traln, branch22states2score insertions, LikelihoodEvaluator& eval) const
{
  auto result = weightMap{}; 
  double minWeight = std::numeric_limits<double>::max(); 

  if(_parallelReduceAtEnd)
    {
      auto data = std::vector<parsimonyNumber>{}; 
      data.reserve(insertions.size()  * 2 );
      for(auto &pair : insertions)
	{
	  auto &values = std::get<1>(pair); 
	  for(auto& v : values)
	    data.push_back(v);
	}

      data = eval.getParallelSetup().getChainComm().allReduce<>(data);

      nat ctr = 0; 
      for(auto &pair :insertions)
	{
	  auto &values = std::get<1>(pair); 
	  for(auto &v : values)
	    v = data[ctr++]; 
	}
      
    }

  for(auto &elem : insertions)
    {
      double score = 0; 

      auto& scoreArray= std::get<1>(elem); 
      for(nat i = 0 ; i < scoreArray.size() ; ++i)
	score += - ( factors[i]  * _parsWarp) * scoreArray[i]; 

      if(score < minWeight)
	minWeight = score; 

      result[std::get<0>(elem)] = score; 
    }

  double sum = 0; 
  for(auto & elem : result)
    {
      double normalizedWeight =exp(minWeight - elem.second); 
      sum += normalizedWeight; 
      elem.second = normalizedWeight; 
    }

  for(auto &elem : result)
    std::get<1>(elem) /= sum; 

  return result; 
} 



BranchPlain ParsimonySPR::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const
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


branch22states2score ParsimonySPR::determineScoresOfInsertions(TreeAln& traln, BranchPlain prunedTree, Randomness &rand, const branch22states2score &alreadyComputed  )
{
  auto result = branch22states2score{}; 

  nodeptr
    p = prunedTree.findNodePtr(traln), 
    pn = p->next->back , 
    pnn = p->next->next->back ; 

  // prune the subtree 
  traln.clipNode( pn, pnn); 
  p->next->back = p->next->next->back = NULL; 

  // fetch all parsimony scores   
  if(not traln.isTipNode(pn)) 
    {
      testInsertParsimony(traln, pn->next->back, p, result, _depth, alreadyComputed);
      testInsertParsimony(traln, pn->next->next->back, p, result, _depth, alreadyComputed); 
    }
  if(not traln.isTipNode(pnn))
    {
      testInsertParsimony(traln, pnn->next->back,p, result, _depth, alreadyComputed); 
      testInsertParsimony(traln, pnn->next->next->back,p, result, _depth, alreadyComputed); 
    }

  traln.clipNode( p->next, pn ); 
  traln.clipNode( p->next->next, pnn); 

  return result; 
} 


void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{ 
  auto blParams = getBranchLengthsParameterView(); 
  auto prunedTree = determinePrimeBranch(traln, rand); 

#ifdef PRINT_DEBUG_PARS 
  tout << "prime branch =" << prunedTree << std::endl; 
#endif

  _pEval.evaluate(traln, prunedTree.findNodePtr(traln), true );

  // decide upon an spr move 
  auto alreadyComputed = branch22states2score {}; 
  alreadyComputed = determineScoresOfInsertions(traln, prunedTree, rand, alreadyComputed); // #, hastings, prior
  auto weightedInsertions = getWeights(traln, alreadyComputed, eval); 

#ifdef PRINT_DEBUG_PARS
  for(auto &v : weightedInsertions)
    tout << "(" << std::get<0>(v) <<  "," << std::get<1>(v) << ")" << std::endl; 
#endif

  auto r = rand.drawRandDouble01(); 
  auto chosen = std::pair<BranchPlain,double>{}; 
  for(auto &v : weightedInsertions)
    {
      if(r < v.second)
	{
	  chosen = v; 
	  break; 
	}
      else 
	r -= v.second; 
    }

#ifdef PRINT_DEBUG_PARS
  tout << "chose " << std::get<0>(chosen) << std::endl;
#endif

  // determine the branch, we pruned from 
  auto desc = traln.getDescendents(prunedTree);
  auto prunedFromBranch = BranchPlain{std::get<0>(desc).getSecNode() , std::get<1>(desc).getSecNode()}; 

#ifdef PRINT_DEBUG_PARS
  tout << " prime branch was pruned from " << prunedFromBranch << " (with neighbors: " << std::get<0>(desc) << " and " << std::get<1>(desc)<< ")"  << std::endl; 
#endif

  // important: save the move 
  _move.extractMoveInfo(traln, std::make_tuple(prunedTree,chosen.first), getSecondaryParameterView() );

  // auto blBackLogProb = 0; 
  // if(_proposePrimeBranch)
  //   {
  //     auto b = _move.getPathHandle().at(1); 
  //     eval.evaluate(traln, b , false );

  //     // TODO 
  //     assert(blParams.size() == 1 ); 

  //     auto blo = BranchLengthOptimizer(traln, b, 30,    eval.getParallelSetup().getChainComm() , blParams); 
  //     blo.optimizeBranches(traln); 
      
  //     auto optParams = blo.getOptimizedParameters(); 
  //     auto optParam = optParams[0]; 

  //     auto gP = optParam.getProposerDistribution<GammaProposer>(traln, 1.61, 2.); // MEH 
  //     auto absLen = traln.getBranch(b, blParams[0]).getInterpretedLength(traln, blParams[0]) ; 
  //     auto blBackLogProb = gP.getLogProbability(absLen); 

  //     // tout << "probability of proposing " << MAX_SCI_PRECISION << absLen << " is " << SOME_FIXED_PRECISION << backProb << " for " << optParam << std::endl; 
  //   }
  
  _move.applyToTree(traln, getSecondaryParameterView() ); 
  
  ParsimonyEvaluator::disorientNode(prunedTree.findNodePtr(traln)); 
  _pEval.evaluateSubtree(traln, prunedTree.findNodePtr(traln));
  
  auto backScores = determineScoresOfInsertions(traln, prunedTree,rand, alreadyComputed); 
  auto weightedInsertionsBack = getWeights(traln, backScores, eval); 
  
  assert(weightedInsertionsBack.find(prunedFromBranch) != weightedInsertionsBack.end()); 
  auto backProb = weightedInsertionsBack[prunedFromBranch]; 
  auto forwProb = std::get<1>(chosen);


  // auto blForwLogProb = 0; 
  // if(_proposePrimeBranch)
  //   {
  //     // now the actual proposing 
      
  //     // tout << SHOW(_move) << std::endl; 

  //     auto invMove = _move.getInverseMove(traln, blParams); 
  //     // tout << SHOW(invMove) << std::endl; 
  //     auto b = invMove.getPathHandle().at(1); 

  //     eval.evaluate(traln, b, false); 
       
  //     auto oldAbsLen = traln.getBranch(b, blParams[0]).getInterpretedLength(traln, blParams[0]);

  //     auto blo = BranchLengthOptimizer(traln,b, 30, eval.getParallelSetup().getChainComm(), blParams); 
  //     blo.optimizeBranches(traln); 
      
  //     auto optParams = blo.getOptimizedParameters();
  //     auto optParam = optParams[0]; 
      
  //     auto gP = optParam.getProposerDistribution<GammaProposer>(traln,1.61, 2. ); // MEH 
      
  //     auto newBranch = gP.proposeBranch(b, traln, blParams[0], rand); 

  //     traln.setBranch(newBranch, blParams[0]); 
      
  //     auto absLen = newBranch.getInterpretedLength(traln, blParams[0]); 

  //     // prior 
  //     auto prNew = blParams[0]->getPrior()->getLogProb( ParameterContent{{absLen}} ); 
  //     auto prOld = blParams[0]->getPrior()->getLogProb( ParameterContent{{oldAbsLen}} ); 
  //     prior.addToRatio( prNew / prOld ); 

  //     auto  blForwLogProb = gP.getLogProbability(absLen); 

  //     // tout << "proposed " << absLen << " with "<< blBackLogProb - blForwLogProb << std::endl; 
  //   }

#ifdef PRINT_DEBUG_PARS
  for(auto &v : weightedInsertionsBack)
    tout << "(" << std::get<0>(v) << "," << std::get<1>(v) << ")" << std::endl; 
  tout << "backProb=" << backProb << " forwProb=" << forwProb << std::endl; 
#endif

  // hastings *= log_double::fromLog((log(backProb) ) - (log(forwProb)));
  hastings *=   ( log_double::fromAbs(backProb)  /   log_double::fromAbs(forwProb)) ; 
} 


void ParsimonySPR::traverse(const TreeAln &traln, nodeptr p, int distance )
{
  std::cout <<  "[" << distance << "] "; 
  for(int i = 0; i < distance; ++i)
    std::cout << " " ; 
  std::cout << p->number << std::endl; 
  if(not traln.isTipNode(p) )
    {
      traverse(traln, p->next->back, distance+1); 
      traverse(traln, p->next->next->back, distance+1);       
    }    
}


void ParsimonySPR::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{  
  auto toEval = _move.getEvalBranch(traln);

  for(auto &elem : _move.getDirtyNodes(traln, false))
    evaluator.markDirty(traln, elem); 

#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL " << toEval << std::endl; 
#endif

  evaluator.evaluate(traln,toEval, false); 
}


void ParsimonySPR::resetState(TreeAln &traln) 
{
  auto params = getBranchLengthsParameterView();
  _move.revertTree(traln, getSecondaryParameterView() ); 
}
 

void ParsimonySPR::autotune() 
{
  // nothing to do 
}
 

AbstractProposal* ParsimonySPR::clone() const
{
  return new ParsimonySPR(*this); 
} 

void ParsimonySPR::printParams(std::ostream &out)  const 
{
  out << " ; radius=" << _depth ; 
}


std::vector<nat> ParsimonySPR::getInvalidatedNodes(const TreeAln& traln) const
{
  return _move.getDirtyNodes(traln, false); 
} 
