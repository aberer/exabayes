#include "ParsimonySPR.hpp"
#include "Topology.hpp"
#include "eval.h"
#include "treeRead.h"

#define DEBUG_PARS_SPR

class  InsertionScore
{
public: 
  InsertionScore(branch _b, vector<nat> _tmp) : b(_b), partitionParsimony(_tmp){}  
  branch getBranch() const  {return b; }

  double getWeight() const {return  logProb; }
  void setWeight(double w) { logProb = w; }

  nat getScore() const
  {
    nat result = 0; 
    for(auto b : partitionParsimony)
      result += b; 
    return result; 
  }

  nat getPartitionScore(int model) const{return partitionParsimony[model] ; }


private: 
  branch b; 
  vector<nat> partitionParsimony; 
  double logProb;
  

  friend ostream& operator<< (ostream &out, const InsertionScore &rhs) { 
    out <<  "(" << rhs.b.thisNode << "," << rhs.b.thatNode << ")" << "=" ; 
    for(auto elem : rhs.partitionParsimony)
      out << elem << "," ; 
    return out; 
  }
} ; 




ParsimonySPR::ParsimonySPR( double _relativeWeight, double _parsWarp, double _blMulti)
  : parsWarp(_parsWarp)    
  , blMulti(_blMulti)    
{
  this->name = "parsSPR"; 
  this->category = TOPOLOGY ; 
  this->relativeProbability =  _relativeWeight; 
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
  // cout << "pushing " << i << endl; 
  
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


void ParsimonySPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
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
  cout << initScore << endl; 

  branch prunedTree  = rand.drawSubtreeUniform(traln);  

  // prune the subtree 
  nodeptr p = findNodeFromBranch(tr, prunedTree),
    pn = p->next->back,
    pnn = p->next->next->back;   
  traln.clipNodeDefault( pn, pnn); 
  p->next->back = p->next->next->back = NULL; 
  branch pruningBranch = constructBranch(pn->number, pnn->number); 

  cout << "pruned  " << p->number << " from "  << pruningBranch << endl;  

  // fetch all parsimony scores   
  if(not traln.isTipNode(pn)) 
    {
      cout << "descending on " << pn->number << " with neighbors " << pn->next->back->number << "," << pn->next->next->back->number  << endl; 
      testInsertParsimony(traln, pn->next->back, p, insertionPoints);
      testInsertParsimony(traln, pn->next->next->back, p, insertionPoints); 
    }
  if(not traln.isTipNode(pnn))
    {
      cout << "descending on " << pnn->number << " with neighbors " << pnn->next->back->number << "," << pnn->next->next->back->number  << endl; 
      testInsertParsimony(traln, pnn->next->back,p, insertionPoints); 
      testInsertParsimony(traln, pnn->next->next->back,p, insertionPoints); 
    }

  // probably not even necessary 
  traln.clipNodeDefault( p->next, pn ); 
  traln.clipNodeDefault( p->next->next, pnn); 
  
  for(auto elem : insertionPoints)
    {
      cout << elem << "," << endl; 
    }

  
  double minWeight = numeric_limits<double>::max(); 
  // now get the real scores 
  for(InsertionScore& elem : insertionPoints)
    {     
      double result = 0; 
      for(int i = 0 ; i < traln.getNumberOfPartitions(); ++i)
	{
	  double states  = double(traln.getPartition(i)->states); 
	  double divFactor = - (parsWarp *  log((1.0/states) - exp(-(states/(states-1) * 0.05 * tr->fracchange)) / states)) ; 
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
      cout << elem.getWeight() << endl; 
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
  cout << "will insert " << p->number << "," << p->back->number << " into " << chosen << endl; 

  path.findPath(traln ,p, findNodeFromBranch(traln.getTr(), chosen.getBranch()));
  path.reverse();  
  
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

  // TODO this is crazy! remove all the redundancies later.
  cout << path << endl; 

  

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
