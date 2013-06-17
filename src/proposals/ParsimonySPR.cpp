#include <unordered_map>
#include <functional>

#include "ParsimonySPR.hpp"
#include "Topology.hpp"
#include "eval.h"
#include "treeRead.h"
#include "output.h"
#include "InsertionScore.hpp"
#include "Branch.hpp"
// #include "PossibilityCollection.hpp"




// #define DEBUG_PARS_SPR


ParsimonySPR::ParsimonySPR(  double _parsWarp, double _blMulti)
  : parsWarp(_parsWarp)    
  , blMulti(_blMulti)    
{
  this->name = "parsSPR"; 
  this->category = TOPOLOGY ; 
  relativeWeight = 5.;

  // cout << "initialized parsiminy spr with warp "<< parsWarp << endl;  
}


static void testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, unordered_map<Branch,vector<nat>, BranchHashNoLength, BranchEqualNoLength > &posses)
{
  nodeptr insertBack =  insertPos->back;   
  traln.clipNodeDefault(insertPos, prunedTree->next);
  traln.clipNodeDefault( insertBack, prunedTree->next->next); 

  Branch b(insertPos->number, insertBack->number); 

  exa_newViewParsimony(traln, prunedTree);
  
  vector<nat> partitionParsimony ; 
  exa_evaluateParsimony(traln, prunedTree->back, FALSE, partitionParsimony);

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


static void verifyParsimony(TreeAln &traln, nodeptr pruned, InsertionScore &score)
{
#ifdef DEBUG_LNL_VERIFY
  globals.debugTree = &traln; 
  cout << endl; 
  cout << "orig " << traln << endl; 
  cout << "copy " << *(globals.debugTree) << endl; 

  tree *tr = globals.debugTree->getTr(); 

  nodeptr p = findNodeFromBranch(tr, constructBranch(pruned->number,pruned->back->number )),
    pn = findNodeFromBranch(tr, score.getBranch()),
    pnn = findNodeFromBranch(tr, invertBranch(score.getBranch())); 
  traln.clipNodeDefault( p->next->back, p->next->next->back); 
  
  traln.clipNodeDefault( pn, p->next); 
  traln.clipNodeDefault( pnn, p->next->next);  

  vector<nat> partitionParsimony; 
  exa_evaluateParsimony(*(globals.debugTree), globals.debugTree->getTr()->nodep[1]->back, TRUE, partitionParsimony);   
  nat verifiedScore = 0; 
  for(auto elem : partitionParsimony)
    verifiedScore += elem; 

  assert(verifiedScore == score.getScore()); 
#endif
}



weightMap ParsimonySPR::getWeights(const TreeAln& traln, const scoreMap &insertions) const
{
  weightMap result; 
  double minWeight = numeric_limits<double>::max(); 

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
  vector<branch> branches; 
  extractBranches(traln, branches); 
  
  tree *tr = traln.getTr(); 

  scoreMap possibilities; 

  nat initScore = 0; 
  vector<nat> partitionParsimony; 
  exa_evaluateParsimony(traln, traln.getTr()->start, TRUE, partitionParsimony);
  for(auto elem : partitionParsimony)
    initScore += elem; 

  // BAD 
  branch prunedTree; 
  nodeptr p = nullptr, pn = nullptr , pnn = nullptr ; 

  do 
    {
      Branch b  = rand.drawBranchWithInnerNode(traln); 
      prunedTree = b.toLegacyBranch(); 

      p = findNodeFromBranch(tr, prunedTree); 
      pn = p->next->back; 
      pnn = p->next->next->back;         

    } while( (traln.isTipNode(pn) &&  traln.isTipNode(pnn))     ); 

  // cout <<"pruning " << p->number << "," << p->back->number << endl; 

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
  pair<Branch,double> chosen; 
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

  path.clear();   
  path.findPath(traln, p, chosen.first.findNodeFromBranch(traln));
  path.reverse();   
  
  assert(possibilities.find(initBranch) == possibilities.end()); 
  possibilities[initBranch] = partitionParsimony; 
  possibilities.erase(chosen.first); 

  auto backWeights = getWeights(traln,possibilities); 
  
  updateHastings( hastings, backWeights[initBranch] /  chosen.second , name); 

  nodeptr aNode = findNodeFromBranch(traln.getTr(), path.at(0)); 
  if(traln.isTipNode(aNode))
    aNode = findNodeFromBranch(traln.getTr(), invertBranch(path.at(0)));
  assert(not traln.isTipNode(aNode)); 
  branch b = path.at(0); 
  while(branchEqualUndirected(b, path.at(0)) || branchEqualUndirected(b,constructBranch(p->number, p->back->number))) 
    {
      aNode = aNode->next; 
      b = constructBranch(aNode->number, aNode->back->number) ;
    }
  path.reverse();
  path.append(b);  
  path.reverse();

#ifdef EFFICIENT 
  // branch mapping would be more efficient
  assert(0); 
#endif

  for(auto &b : branches)    
    {
      Branch bN ; 
      bN.initFromLegacy(b); 
      auto p = bN.findNodeFromBranch(traln); 
      traln.clipNode(p,p->back, b.length[0]); 
    }
} 


void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{ 
  path.clear(); 
  determineSprPath(traln, rand, hastings, prior); 
  path.saveBranchLengthsPath(traln);  

  move.applyPathAsESPR(traln,path);

  bool modifiesBl = false; 
  for(auto v : secVar)
    modifiesBl |= v.getCategory() == BRANCH_LENGTHS; 

  assert(traln.getNumBranches() == 1); 


#ifdef  NO_SEC_BL_MULTI
  modifiesBl = false;
#endif

  if( modifiesBl)
    {
      auto brPr =  secVar[0].getPrior();
      move.multiplyAlongBranchESPR(traln, rand, hastings, prior, path, blMulti, brPr);      
    }

  debug_checkTreeConsistency(traln); 
}



static void traverse(const TreeAln &traln, nodeptr p, int distance )
{
  cout <<  "[" << distance << "] "; 
  for(int i = 0; i < distance; ++i)
    cout << " " ; 
  cout << p->number << endl; 
  if(not traln.isTipNode(p) )
    {
      traverse(traln, p->next->back, distance+1); 
      traverse(traln, p->next->next->back, distance+1);       
    }    
}


void ParsimonySPR::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{  
  tree *tr = traln.getTr(); 

  branch b = getThirdBranch(tr, path.at(0), path.at(1)); 
  branch lastBranch = path.at(path.size()-1); 

  branch bA = constructBranch(b.thisNode, lastBranch.thisNode), 
    bB = constructBranch(b.thisNode,lastBranch.thatNode); 

  assert(branchExists(tr, bA) && branchExists(tr,bB)); 

  branch futureRoot = getThirdBranch(tr, bA, bB ); 

  nodeptr toEval = findNodeFromBranch(tr, futureRoot);

  move.destroyOrientationAlongPath(path, traln, toEval);

  evaluateGenericWrapper(traln, toEval, FALSE);
}

void ParsimonySPR::resetState(TreeAln &traln, PriorBelief &prior) 
{
  // cout << "RESET" << endl; 
  move.resetAlongPathForESPR (traln, prior, path);   
  path.restoreBranchLengthsPath(traln, prior ); 
  debug_checkTreeConsistency(traln);
  debug_printTree(traln); 
}
 

void ParsimonySPR::autotune() 
{
  // nothing to do 
}
 

AbstractProposal* ParsimonySPR::clone() const
{
  return new ParsimonySPR(*this); 
} 

