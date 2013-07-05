#include <unordered_map>
#include <functional>

#include "ParsimonySPR.hpp"
#include "Topology.hpp"
#include "treeRead.h"
#include "InsertionScore.hpp"
#include "Branch.hpp"


ParsimonySPR::ParsimonySPR(  double _parsWarp, double _blMulti)
  : parsWarp(_parsWarp)    
  , blMulti(_blMulti)    
{
  this->name = "parsSPR"; 
  this->category = Category::TOPOLOGY ; 
  relativeWeight = 5.;
}


void ParsimonySPR::testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, std::unordered_map<Branch,std::vector<nat>, BranchHashNoLength, BranchEqualNoLength > &posses)
{
  nodeptr insertBack =  insertPos->back;   
  traln.clipNodeDefault(insertPos, prunedTree->next);
  traln.clipNodeDefault( insertBack, prunedTree->next->next); 

  Branch b(insertPos->number, insertBack->number); 

  std::vector<nat> partitionParsimony ; 
  pEval.evaluateSubtree(traln,  prunedTree); 
  pEval.evaluateSubtree(traln,  prunedTree->back); // necessary? 
  pEval.evaluate(traln, prunedTree, false, partitionParsimony); 

  assert(posses.find(b) == posses.end()) ;   
  posses[b] = partitionParsimony; 

  traln.clipNodeDefault(insertPos, insertBack); 
  prunedTree->next->back = prunedTree->next->next->back = NULL; 

  // recursively descend 
  if(not traln.isTipNode(insertPos))
    {
      testInsertParsimony(traln, insertPos->next->back, prunedTree, posses); 
      testInsertParsimony(traln, insertPos->next->next->back, prunedTree, posses); 
    }
}


// void ParsimonySPR::verifyParsimony(TreeAln &traln, nodeptr pruned)
// {
// #ifdef DEBUG_LNL_VERIFY
//   globals.debugTree = &traln; 
//   cout << endl; 
//   cout << "orig " << traln << endl; 
//   cout << "copy " << *(globals.debugTree) << endl; 

//   tree *tr = globals.debugTree->getTr(); 

//   nodeptr p = findNodeFromBranch(tr, constructBranch(pruned->number,pruned->back->number )),
//     pn = findNodeFromBranch(tr, score.getBranch()),
//     pnn = findNodeFromBranch(tr, invertBranch(score.getBranch())); 
//   traln.clipNodeDefault( p->next->back, p->next->next->back); 
  
//   traln.clipNodeDefault( pn, p->next); 
//   traln.clipNodeDefault( pnn, p->next->next);  

//   vector<nat> partitionParsimony; 
//   exa_evaluateParsimony(*(globals.debugTree), globals.debugTree->getTr()->nodep[1]->back, TRUE, partitionParsimony);   
//   nat verifiedScore = 0; 
//   for(auto elem : partitionParsimony)
//     verifiedScore += elem; 

//   assert(verifiedScore == score.getScore()); 
// #endif
// }



weightMap ParsimonySPR::getWeights(const TreeAln& traln, const scoreMap &insertions) const
{
  weightMap result; 
  double minWeight = std::numeric_limits<double>::max(); 

  for(auto &elem : insertions)
    {
      double score = 0; 
      for(int i = 0 ; i < traln.getNumberOfPartitions(); ++i)
	{
	  double states  = double(traln.getPartition(i)->states); 
	  double divFactor = - (parsWarp *  log((1.0/states) - exp(-(states/(states-1) * 0.05)) / states)) ;  //  * tr->fracchange
	  score += divFactor * elem.second[i];
	}

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



void ParsimonySPR::determineSprPath(TreeAln& traln, Randomness &rand, double &hastings, PriorBelief &prior )
{
  std::vector<Branch> branches = traln.extractBranches(); 

  scoreMap possibilities; 

  std::vector<nat> partitionParsimony; 
  pEval.evaluate(traln, traln.getTr()->start, true, partitionParsimony);

  // BAD 
  Branch prunedTree; 
  nodeptr p = nullptr, pn = nullptr , pnn = nullptr ; 

  do 
    {
      prunedTree  = TreeRandomizer::drawBranchWithInnerNode(traln,rand); 
      p = prunedTree.findNodePtr(traln);
      pn = p->next->back; 
      pnn = p->next->next->back;         

    } while( (traln.isTipNode(pn) &&  traln.isTipNode(pnn))     ); 

  Branch initBranch = Branch(pn->number, pnn->number); 
  
  // prune the subtree 
  traln.clipNodeDefault( pn, pnn); 
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
  traln.clipNodeDefault( p->next, pn ); 
  traln.clipNodeDefault( p->next->next, pnn); 

  auto weightedInsertions = getWeights(traln, possibilities) ; 

  double r = rand.drawRandDouble01(); 
  std::pair<Branch,double> chosen; 
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


  {
#ifdef EFFICIENT 
    // branch mapping would be more efficient
    assert(0); 
#endif
    for(auto &b : branches)    
      {
	auto p = b.findNodePtr(traln); 
	double tmp = b.getLength(); 
	traln.clipNode(p,p->back, tmp); 
      }
  }

  // important: save the move 
  Branch pruned(p->number,p->back->number); 

  move.extractMoveInfo(traln, {pruned,chosen.first});

  // update the hastings  
  assert(possibilities.find(initBranch) == possibilities.end()); 
  possibilities[initBranch] = partitionParsimony; 
  possibilities.erase(chosen.first); 
  auto backWeights = getWeights(traln,possibilities);   
  updateHastings( hastings, backWeights[initBranch] /  chosen.second , name); 
} 


void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{ 
  determineSprPath(traln, rand, hastings, prior); 
  move.applyToTree(traln); 

  bool modifiesBl = false; 
  for(auto &v : secVar)
    modifiesBl |= v->getCategory() == Category::BRANCH_LENGTHS; 

  assert(traln.getNumBranches() == 1); 

#ifdef  NO_SEC_BL_MULTI
  modifiesBl = false;  
  // assert(0);
#endif

  if( modifiesBl)
    {
      // assert(0); 
      auto brPr =  secVar[0]->getPrior();
      move.multiplyBranches(traln, rand, hastings, prior,  blMulti,{ brPr}); 
    }

  // debug_checkTreeConsistency(traln);   
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


void ParsimonySPR::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior) 
{  
  Branch toEval = move.getEvalBranch(traln);
  nodeptr toEvalP = toEval.findNodePtr(traln) ; 
  move.disorientAtNode(traln,toEvalP);
  evaluator.evaluate(traln,toEval, false); 
}

void ParsimonySPR::resetState(TreeAln &traln, PriorBelief &prior) 
{
  move.revertTree(traln,prior); 
}
 

void ParsimonySPR::autotune() 
{
  // nothing to do 
}
 

AbstractProposal* ParsimonySPR::clone() const
{
  return new ParsimonySPR(*this); 
} 
