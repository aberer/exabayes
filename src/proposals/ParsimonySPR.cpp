#include "ParsimonySPR.hpp"
#include "Topology.hpp"
#include "eval.h"
#include "treeRead.h"
#include "output.h"
#include "InsertionScore.hpp"
#include "Branch.hpp"


// #define DEBUG_PARS_SPR


ParsimonySPR::ParsimonySPR(  double _parsWarp, double _blMulti)
  : parsWarp(_parsWarp)    
  , blMulti(_blMulti)    
{
  this->name = "parsSPR"; 
  this->category = TOPOLOGY ; 
  relativeWeight = 5.;
}


static void testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, vector<InsertionScore> &insertPoints)
{
  nodeptr insertBack =  insertPos->back;   
  traln.clipNodeDefault(insertPos, prunedTree->next);
  traln.clipNodeDefault( insertBack, prunedTree->next->next); 
  branch b = constructBranch(insertPos->number, insertBack->number); 

  exa_newViewParsimony(traln, prunedTree);
  
  vector<nat> partitionParsimony ; 
  exa_evaluateParsimony(traln, prunedTree->back, FALSE, partitionParsimony);

  InsertionScore i(b, partitionParsimony); 
  insertPoints.push_back(i); 

  traln.clipNodeDefault(insertPos, insertBack); 
  prunedTree->next->back = prunedTree->next->next->back = NULL; 

  // recursively descend 
  if(not traln.isTipNode(insertPos))
    {
      testInsertParsimony(traln, insertPos->next->back, prunedTree, insertPoints); 
      testInsertParsimony(traln, insertPos->next->next->back, prunedTree, insertPoints); 
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





void ParsimonySPR::determineSprPath(TreeAln& traln, Randomness &rand, double &hastings, PriorBelief &prior )
{
  vector<branch> branches; 
  extractBranches(traln, branches); 
  
  vector<InsertionScore> insertionPoints; 
  tree *tr = traln.getTr(); 

  nat initScore = 0; 
  vector<nat> partitionParsimony; 
  exa_evaluateParsimony(traln, traln.getTr()->start, TRUE, partitionParsimony);
  for(auto elem : partitionParsimony)
    initScore += elem; 

  // BAD 
  branch prunedTree; 
  nodeptr p = nullptr, pn = nullptr , pnn = nullptr ; 

  while( ( pn == nullptr && pnn == nullptr ) 
	 || ( traln.isTipNode(pn) && traln.isTipNode(pnn) ) )
    {      
      prunedTree  = rand.drawSubtreeUniform(traln);  
      p = findNodeFromBranch(tr, prunedTree); 
      pn = p->next->back; 
      pnn = p->next->next->back;   
    }

  branch initBranch = constructBranch(pn->number, pnn->number); 
  insertionPoints.push_back(InsertionScore(initBranch, partitionParsimony)); 

  // prune the subtree 
  traln.clipNodeDefault( pn, pnn); 
  p->next->back = p->next->next->back = NULL; 

  // fetch all parsimony scores   
  if(not traln.isTipNode(pn)) 
    {
      testInsertParsimony(traln, pn->next->back, p, insertionPoints);
      testInsertParsimony(traln, pn->next->next->back, p, insertionPoints); 
    }
  if(not traln.isTipNode(pnn))
    {
      testInsertParsimony(traln, pnn->next->back,p, insertionPoints); 
      testInsertParsimony(traln, pnn->next->next->back,p, insertionPoints); 
    }

  // probably not even necessary 
  traln.clipNodeDefault( p->next, pn ); 
  traln.clipNodeDefault( p->next->next, pnn); 


  double minWeight = numeric_limits<double>::max(); 
  // now get the real scores 
  for(InsertionScore& elem : insertionPoints)
    {     
      double result = 0; 
      for(int i = 0 ; i < traln.getNumberOfPartitions(); ++i)
	{
	  double states  = double(traln.getPartition(i)->states); 
	  double divFactor = - (parsWarp *  log((1.0/states) - exp(-(states/(states-1) * 0.05)) / states)) ;  //  * tr->fracchange
	  result += divFactor * elem.getPartitionScore(i);
	}

      if(result < minWeight)
	minWeight = result; 

      elem.setWeight(result); 
    }

  double sum =  0; 
  for(InsertionScore& elem : insertionPoints)
    {
      double normWeight = exp( minWeight - elem.getWeight())  ; 
      sum += normWeight;       
      elem.setWeight(normWeight); 
    }
  
  // loop could be avoided 
  for(InsertionScore & elem : insertionPoints)
    {
      double val = elem.getWeight(); 
      elem.setWeight(val / sum); 
    }

  // draw the proposal 
  double r = rand.drawRandDouble01();
  InsertionScore& chosen = insertionPoints[0];

  Branch initBranchNew; 
  initBranchNew.initFromLegacy(initBranch); 

  do 
    {
      for(InsertionScore &elem : insertionPoints)
	{
	  double weight = elem.getWeight(); 
	  if(r < weight)
	    {
	      chosen = elem; 
	      break; 
	    }
	  else 
	    r -= weight; 
	}
    }while (chosen.getNewBranch().equalsUndirected(initBranchNew) ); 

  path.clear();   
  path.findPath(traln ,p, findNodeFromBranch(traln.getTr(), chosen.getBranch()) );
  path.reverse(); 


  // adjust hastings   
  double initWeight = 0; 
  for(auto &elem: insertionPoints)
    {
      if(elem.getNewBranch().equalsUndirected(initBranchNew))
	initWeight = elem.getWeight(); 
    }
  assert(initWeight != 0); 

  updateHastings(hastings, initWeight / chosen.getWeight(), name); 

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

