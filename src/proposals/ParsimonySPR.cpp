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
  
  cout << "trying to find " << constructBranch(pruned->number,pruned->back->number ) << ", "
       << score.getBranch() <<  ","
       << invertBranch(score.getBranch()) << endl; 
  
  nodeptr p = findNodeFromBranch(tr, constructBranch(pruned->number,pruned->back->number )),
    pn = findNodeFromBranch(tr, score.getBranch()),
    pnn = findNodeFromBranch(tr, invertBranch(score.getBranch())); 
  traln.clipNodeDefault( p->next->back, p->next->next->back); 
  
  traln.clipNodeDefault( pn, p->next); 
  traln.clipNodeDefault( pnn, p->next->next);  

  cout << "after spr " << *(globals.debugTree) << endl ; 
  
  vector<nat> partitionParsimony; 
  exa_evaluateParsimony(*(globals.debugTree), globals.debugTree->getTr()->nodep[1]->back, TRUE, partitionParsimony);   
  nat verifiedScore = 0; 
  for(auto elem : partitionParsimony)
    verifiedScore += elem; 
  
  // cout << "for " << score.getBranch() << "\t toVerify: " << score.getScore() <<  "\tveification result: " << verifiedScore <<  endl; 
  assert(verifiedScore == score.getScore()); 
#endif
}





void ParsimonySPR::determineSprPath(TreeAln& traln, Randomness &rand, double &hastings, PriorBelief &prior )
{
  vector<InsertionScore> insertionPoints; 
  tree *tr = traln.getTr(); 
  Topology prevTopo(traln.getTr()->mxtips) ; 
  prevTopo.saveTopology(traln); 

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


#ifdef UNSURE
  // dont choose the original position again 
  assert(0); 
#endif  

  // draw the proposal 
  double r = rand.drawRandDouble01();
  InsertionScore& chosen = insertionPoints[0]; 
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

  path.clear();   

  Branch chosenBranch;   
  path.findPath(traln ,p, findNodeFromBranch(traln.getTr(), chosen.getBranch()) );

  cout << "now path is " << path << endl; 

  cout << "will insert " << p->number << "," << p->back->number 
       << " currently hooked to " << p->next->back->number << "," << p->next->next->back->number
       << " into " << chosen << endl; 
  
  cout << "current root is " << findRoot(traln.getTr() ) << endl; 

  int numTax = traln.getTr()->mxtips; 

  for(int i = 1 ; i < numTax; ++i)
    assert(traln.getTr()->nodep[i]->x == 0 ); 

  path.reverse();  

  // TODO : here's the bug!   
  
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

  cout << "path is " << path << endl; 

  // TODO this is crazy! remove all the redundancies later.
  // cout << path << endl; 
  
 
#ifdef UNSURE
  // modify the hastings...not today
  assert(0); 
#endif 
  
  prevTopo.restoreTopology(traln);
} 


void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{ 

  path.clear(); 
  determineSprPath(traln, rand, hastings, prior); 
  path.saveBranchLengthsPath(traln);  

  move.applyPathAsESPR(traln,path);


#ifdef ESPR_MULTIPLY_BL
  bool modifiesBl = false; 
  for(auto v : secVar)
    modifiesBl |= v.getCategory() == BRANCH_LENGTHS; 

  assert(traln.getNumBranches() == 1); 

  if( modifiesBl)
    {
      auto brPr =  secVar[0].getPrior();
      move.multiplyAlongBranchESPR(traln, rand, hastings, prior, path, blMulti, brPr);      
    }
#ifdef DEBUG_SHOW_TOPO_CHANGES
  cout << "after multiply: " << traln  << endl; 
#endif
#endif
  
  debug_checkTreeConsistency(traln); 

  
  TreePrinter tp(false, true, false); 
  cout << "after psr-application " << tp.printTree(traln) << endl; 
  
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

#if 1 
  evaluateGenericWrapper(traln, toEval, FALSE);
  cout << "tree aln: " << traln.getTr()->likelihood << endl; 
#else 
  evaluateGenericWrapper(traln, traln.getTr()->start, TRUE);
#endif
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

