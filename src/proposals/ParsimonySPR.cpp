#include <unordered_map>
#include <functional>

#include "priors/AbstractPrior.hpp"

#include "ParsimonySPR.hpp"
#include "Branch.hpp"

#include "TreePrinter.hpp"


ParsimonySPR::ParsimonySPR( double parsWarp, double blMulti)
  : AbstractProposal(Category::TOPOLOGY, "parsSPR", 5. , false)
  , _parsWarp(parsWarp)    
  , _blMulti(blMulti)    
{
}


void ParsimonySPR::testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, std::unordered_map<BranchPlain,std::vector<nat> > &posses)
{
  nodeptr insertBack =  insertPos->back;   
  traln.clipNode(insertPos, prunedTree->next);
  traln.clipNode( insertBack, prunedTree->next->next); 

  auto b = BranchPlain(insertPos->number, insertBack->number); 

  auto partitionParsimony = std::vector<nat> {}; 
  auto pLengthAtBranch = std::vector<nat> {}; 
  _pEval.evaluateSubtree(traln, prunedTree); 
  _pEval.evaluateSubtree(traln, prunedTree->back); // necessary? 
  _pEval.evaluate(traln, prunedTree, false, partitionParsimony, pLengthAtBranch); 

  assert(posses.find(b) == posses.end()) ;   
  posses[b] = partitionParsimony; 

  traln.clipNode(insertPos, insertBack); 
  prunedTree->next->back = prunedTree->next->next->back = NULL; 

  // recursively descend 
  if(not traln.isTipNode(insertPos))
    {
      testInsertParsimony(traln, insertPos->next->back, prunedTree, posses); 
      testInsertParsimony(traln, insertPos->next->next->back, prunedTree, posses); 
    }
}


weightMap ParsimonySPR::getWeights(const TreeAln& traln, const scoreMap &insertions) const
{
  auto result = weightMap{}; 
  double minWeight = std::numeric_limits<double>::max(); 

  auto factors = std::vector<double>(traln.getNumberOfPartitions()); 
  double typicalBL = 0.05; 
  for(nat i = 0; i < traln.getNumberOfPartitions( ); ++i )
    {
      auto& partition = traln.getPartition(i); 
      nat states = partition.states; 
      factors[i] =  - (_parsWarp * log((1.0/states) - exp(-(states/(states-1) * typicalBL)) / states));
    }	  
  
  for(auto &elem : insertions)
    {
      double score = 0; 
      for(nat i = 0 ; i < traln.getNumberOfPartitions(); ++i)
	score += factors[i] * elem.second[i];

      if(score < minWeight)
	minWeight = score; 

      result[elem.first] = score; 
    }

  double sum = 0; 
  for(auto & elem : result)
    {
      double normalizedWeight =exp(minWeight - elem.second); 
      sum += normalizedWeight; 
      elem.second = normalizedWeight; 
    }

  for(auto &elem : result)
    elem.second /= sum; 

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


void ParsimonySPR::determineSprPath(TreeAln& traln, Randomness &rand, double &hastings, PriorBelief &prior )
{
  auto branches = traln.extractBranches( ); 
  auto possibilities = scoreMap{}; 

  auto partitionParsimony =  std::vector<nat>{}; 
  auto pLengthAtBranch =  std::vector<nat>{}; 
  _pEval.evaluate(traln, traln.getAnyBranch().findNodePtr(traln), true, partitionParsimony, pLengthAtBranch);

  auto prunedTree = determinePrimeBranch(traln, rand); 

  // needed? 
  nodeptr
    p = prunedTree.findNodePtr(traln), 
    pn = p->next->back , 
    pnn = p->next->next->back ; 


  auto initBranch = BranchPlain(pn->number, pnn->number); 
  
  // prune the subtree 
  traln.clipNode( pn, pnn); 
  p->next->back = p->next->next->back = NULL; 

  // fetch all parsimony scores   
  if(not traln.isTipNode(pn)) 
    {
      testInsertParsimony(traln, pn->next->back, p, possibilities);
      testInsertParsimony(traln, pn->next->next->back, p, possibilities); 
    }
  if(not traln.isTipNode(pnn))
    {
      testInsertParsimony(traln, pnn->next->back,p, possibilities); 
      testInsertParsimony(traln, pnn->next->next->back,p, possibilities); 
    }

  // probably not even necessary 
  traln.clipNode( p->next, pn ); 
  traln.clipNode( p->next->next, pnn); 

  auto weightedInsertions = getWeights(traln, possibilities) ; 

  double r = rand.drawRandDouble01(); 
  std::pair<BranchPlain,double> chosen; 
  for(auto v : weightedInsertions)
    {
      if(r < v.second)
	{
	  chosen = v; 
	  break; 
	}
      else 
	r -= v.second; 
    }

  // important: save the move 
  auto pruned = BranchPlain(p->number,p->back->number); 

  _move.extractMoveInfo(traln, std::make_tuple(pruned,chosen.first), getSecondaryParameterView() );

  // update the hastings  
  assert(possibilities.find(initBranch) == possibilities.end()); 
  possibilities[initBranch] = partitionParsimony; 
  possibilities.erase(chosen.first); 
  auto backWeights = getWeights(traln,possibilities);   
  AbstractProposal::updateHastingsLog( hastings, log(backWeights[initBranch]) - log(chosen.second) , _name); 
} 



void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{ 
  auto blParams = getBranchLengthsParameterView(); 

  determineSprPath(traln, rand, hastings, prior); 

  _move.applyToTree(traln, getSecondaryParameterView() ); 

  bool modifiesBl = false; 
  for(auto &v : _secondaryParameters)
    modifiesBl |= v->getCategory() == Category::BRANCH_LENGTHS; 

#ifdef  NO_SEC_BL_MULTI
  modifiesBl = false;  
  // assert(0);
#endif

  if( modifiesBl)
    {
      assert(0); 
    }
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
  // auto tp =  TreePrinter{false, true, false};

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


std::vector<nat> ParsimonySPR::getInvalidatedNodes(const TreeAln& traln) const
{
  return _move.getDirtyNodes(traln, false); 
} 
