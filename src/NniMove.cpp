#include "NniMove.hpp"
#include "LikelihoodEvaluator.hpp"
#include "AbstractProposal.hpp"


void NniMove::init(const TreeAln &traln, const Branch& _innerBranch, pair<int,int> _switching)
{
  outerBranches.clear(); 

  innerBranch = _innerBranch; 
  switching = _switching; 
  nodeptr p = innerBranch.findNodeFromBranch(traln); 
  assert(traln.getNumBranches() == 1 ); 
  innerBranch.setLength(p->z[0]); 
  
  outerBranches.push_back(Branch(p->next->back->number, p->number, p->z[0])); 
  outerBranches.push_back(Branch(p->next->next->back->number, p->number, p->z[0])); 
  outerBranches.push_back(Branch(p->back->next->back->number,p->back->number, p->back->z[0])); 
  outerBranches.push_back(Branch(p->back->next->next->back->number, p->back->number, p->back->z[0])); 

  bool foundA = false; 
  bool foundB = false; 
  for(auto &v : outerBranches)
    {
      foundA |= v.nodeIsInBranch(switching.first); 
      foundB |= v.nodeIsInBranch(switching.second); 
    }
}


void NniMove::disortient(TreeAln &traln) const 
{
  nodeptr p = innerBranch.findNodeFromBranch(traln); 
  LikelihoodEvaluator::disorientNode(p); 
  LikelihoodEvaluator::disorientNode(p->back); 
}



void NniMove::apply(TreeAln &traln) const 
{  
  Branch a(switching.first, innerBranch.getPrimNode()); 
  Branch b; 
  if( not a.exists(traln))
    {
      a = Branch(switching.first, innerBranch.getSecNode()); 
      assert(a.exists(traln)); 
      b = Branch(switching.second, innerBranch.getPrimNode()); 
      assert(b.exists(traln));       
    }
  else 
    {
      b = Branch(switching.second, innerBranch.getSecNode()); 
      assert(b.exists(traln));       
    }  

  nodeptr pA = a.findNodeFromBranch(traln),
    pB = b.findNodeFromBranch(traln), 
    pABack = pA->back,
    pBBack = pB->back; 

  traln.clipNode(pA,pBBack, pA->z[0]); 
  traln.clipNode(pB, pABack, pB->z[0]);
}


void NniMove::revert(TreeAln &traln) const
{
  apply(traln);

  // revert branch lengths 
  nodeptr p = innerBranch.findNodeFromBranch(traln); 
  double tmp = innerBranch.getLength(); 
  traln.clipNode(p,p->back, tmp); 

  for(auto &elem : outerBranches)
    {
      nodeptr q = elem.findNodeFromBranch(traln); 
      tmp =  elem.getLength(); 
      traln.clipNode(q,q->back, tmp); 
    }
} 


void NniMove::multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter , vector<shared_ptr<AbstractPrior> > priors, string name)
{
  nodeptr p = branch.findNodeFromBranch(traln); 

  double multi = rand.drawMultiplier(parameter); 
  double oldV = traln.getBranchLength( p,0); 
  double newV = pow(oldV, multi); 

  traln.setBranchLengthBounded(newV, 0, p); 
  
  double realM = log(newV) /  log(oldV);   
  updateHastings(hastings, realM, name); 
  
  assert(priors.size() == 1); 
  auto brPr = priors[0]; 
  prior.updateBranchLengthPrior(traln,oldV, newV, brPr);
}

void NniMove::multiplyAllBranches( TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter, vector<shared_ptr<AbstractPrior> > priors, string name) 
{
  multiplyBranch(innerBranch, traln, hastings, prior, rand, parameter,priors, name); 
  for(auto &v : outerBranches)
    multiplyBranch(v,traln, hastings, prior, rand, parameter, priors, name);     
} 
