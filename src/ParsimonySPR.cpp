#include "ParsimonySPR.hpp"
#include "Topology.hpp"
#include "eval.h"
#include "adapters.h"

#define DEBUG_PARS_SPR

class  InsertionScore
{
public: 
  InsertionScore(branch _b, nat _score) : b(_b), score(_score){}
  branch getBranch() {return b; }
  nat getScore(){return score; }

private: 
  branch b; 
  nat score ;   
  
  friend ostream& operator<< (ostream &out, const InsertionScore &rhs) { return out << "(" << rhs.b.thisNode << "," << rhs.b.thatNode << ")" << "="  << rhs.score ;  }
} ; 




ParsimonySPR::ParsimonySPR( double _relativeWeight, double _parsWarp, double _blMulti)
  : parsWarp(_parsWarp)    
  , blMulti(_blMulti)    
{
  this->name = "parsSPR"; 
  this->category = TOPOLOGY ; 
  ptype = PARSIMONY_SPR; 
  this->relativeProbability =  _relativeWeight; 
}


static void testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, vector<InsertionScore> &insertPoints)
{
  // cout << "inserting" << prunedTree->number << "," << prunedTree->back->number <<   " at pos " << insertPos->number << "," <<  insertPos->back->number << endl; 

  tree *tr = traln.getTr(); 

  nodeptr insertBack =  insertPos->back;   
  exa_hookupDefault(tr, insertPos, prunedTree->next);
  exa_hookupDefault(tr, insertBack, prunedTree->next->next); 
  branch b = constructBranch(insertPos->number, insertBack->number); 

  cout << "newview on " << prunedTree->number << endl; 
  exa_newViewParsimony(traln, prunedTree);
  nat result = exa_evaluateParsimony(traln, prunedTree->back, FALSE);
  nat result2 = exa_evaluateParsimony(traln, prunedTree->back, TRUE);
   
  cout << "result1=" << result << ",result2=" << result2 << endl; 

  InsertionScore i(b, result); 
  insertPoints.push_back(i); 
  // cout << "pushing " << i << endl; 
  
  exa_hookupDefault(tr, insertPos, insertBack); 
  prunedTree->next->back = prunedTree->next->next->back = NULL; 

  // recursively descend 
  if(NOT traln.isTipNode(insertPos))
    {
      testInsertParsimony(traln, insertPos->next->back, prunedTree, insertPoints); 
      testInsertParsimony(traln, insertPos->next->next->back, prunedTree, insertPoints); 
    }
}


static void verifyParsimony(TreeAln &traln, nodeptr pruned, InsertionScore &score)
{
#ifdef DEBUG_LNL_VERIFY
  *globals.debugTree = traln; 
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
  exa_hookupDefault(tr, p->next->back, p->next->next->back); 
  
  exa_hookupDefault(tr, pn, p->next); 
  exa_hookupDefault(tr, pnn, p->next->next);  

  cout << "after spr " << *(globals.debugTree) << endl ; 
  
  nat verifiedScore = exa_evaluateParsimony(*(globals.debugTree), globals.debugTree->getTr()->nodep[1]->back, TRUE); 

  cout << "for " << score.getBranch() << "\t toVerify: " << score.getScore() <<  "\tveification result: " << verifiedScore <<  endl; 
  assert(verifiedScore == score.getScore()); 
#endif
}




void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{    
  vector<InsertionScore> insertionPoints; 
  tree *tr = traln.getTr(); 
  Topology prevTopo(traln.getTr()->mxtips) ; 
  prevTopo.saveTopology(traln); 

  double initScore =   exa_evaluateParsimony(traln, traln.getTr()->start, TRUE);
  cout << initScore << endl; 

  branch prunedTree  = rand.drawSubtreeUniform(traln);  
  // path.append(prunedTree); 

  // prune the subtree 
  nodeptr p = findNodeFromBranch(tr, prunedTree),
    pn = p->next->back,
    pnn = p->next->next->back;   
  exa_hookupDefault(tr, pn, pnn); 
  p->next->back = p->next->next->back = NULL; 
  branch pruningBranch = constructBranch(pn->number, pnn->number); 
  // insertionPoints.push_back(InsertionScore(pruningBranch, initScore)); 

  cout << "pruned  " << p->number << " from "  << pruningBranch << endl;  

  // fetch all parsimony scores   
  if(NOT traln.isTipNode(pn)) 
    {
      cout << "descending on " << pn->number << " with neighbors " << pn->next->back->number << "," << pn->next->next->back->number  << endl; 
      testInsertParsimony(traln, pn->next->back, p, insertionPoints);
      testInsertParsimony(traln, pn->next->next->back, p, insertionPoints); 
    }
  if(NOT traln.isTipNode(pnn))
    {
      cout << "descending on " << pnn->number << " with neighbors " << pnn->next->back->number << "," << pnn->next->next->back->number  << endl; 
      testInsertParsimony(traln, pnn->next->back,p, insertionPoints); 
      testInsertParsimony(traln, pnn->next->next->back,p, insertionPoints); 
    }

  // probably not even necessary 
  exa_hookupDefault( tr, p->next, pn ); 
  exa_hookupDefault( tr, p->next->next, pnn); 
  
  for(auto elem : insertionPoints )
    {
      cout << elem << "," ; 
// #ifdef DEBUG_PARS_SPR
//       verifyParsimony(traln,p, elem); 
// #endif
    }
  cout << endl; 

  
  assert(0);   
  prevTopo.restoreTopology(traln);
}


void ParsimonySPR::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  assert(0); 
}

void ParsimonySPR::resetState(TreeAln &traln, PriorBelief &prior) 
{
  assert(0); 
}
 

void ParsimonySPR::autotune() 
{
  // nothing to do 
}


AbstractProposal* ParsimonySPR::clone() const
{
  return new ParsimonySPR( relativeProbability, parsWarp, blMulti);
} 
