#include "NniMove.hpp"
#include "LikelihoodEvaluator.hpp"
#include "AbstractProposal.hpp"

void NniMove::applyToTree(TreeAln &traln) const 
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

  nodeptr pA = a.findNodePtr(traln),
    pB = b.findNodePtr(traln), 
    pABack = pA->back,
    pBBack = pB->back; 

  traln.clipNode(pA,pBBack, pA->z[0]); 
  traln.clipNode(pB, pABack, pB->z[0]);
} 

void NniMove::revertTree(TreeAln &traln, PriorBelief &prior) const 
{
  applyToTree(traln);

  // revert branch lengths 
  nodeptr p = innerBranch.findNodePtr(traln); 
  double tmp = innerBranch.getLength(); 
  traln.clipNode(p,p->back, tmp); 

  for(auto &elem : outerBranches)
    {
      nodeptr q = elem.findNodePtr(traln); 
      tmp =  elem.getLength(); 
      traln.clipNode(q,q->back, tmp); 
    }
} 


/** 
    @brief disorients the nodes for evaluation at node p. 
 */ 
void NniMove::disorientAtNode(TreeAln &traln, nodeptr p) const 
{ 
  if( Branch(p->number, p->back->number).equalsUndirected(innerBranch))
    {      
      disorientHelper(traln, p); 
      disorientHelper(traln, p->back); 
    }
  else if(  innerBranch.nodeIsInBranch(p->number) )
    {
      assert(not innerBranch.nodeIsInBranch( p->back->number)); 
      disorientAtNode(traln, p);
      int node = innerBranch.getOtherNode(p->number); 
      disorientAtNode(traln, Branch(node, p->number).findNodePtr(traln)); 
    } 
  else 
    {
      // disallowing evaluation on any more peripheral pointers. I guess this would not be needed.   
      assert(0); 
    }
} 



void NniMove::extractMoveInfo(const TreeAln &traln, std::vector<Branch> description) 
{
  assert(description.size() == 2 ) ; 
  std::pair<int,int>  _switching = { description.at(0).findNodePtr(traln)->next->back->number , description.at(1).findNodePtr(traln)->number}; 
  Branch _innerBranch = description.at(0); 

  outerBranches.clear(); 

  innerBranch = _innerBranch; 
  switching = _switching; 
  nodeptr p = innerBranch.findNodePtr(traln); 
  assert(traln.getNumBranches() == 1 ); 
  innerBranch.setLength(p->z[0]); 
  
  nodeptr q =p->next->back; 
  outerBranches.push_back(Branch(q->number, q->back->number, q->z[0])); 
  q = p->next->next->back; 
  outerBranches.push_back(Branch(q->number, q->back->number, q->z[0])); 
  q = p->back->next->back; 
  outerBranches.push_back(Branch(q->number, q->back->number, q->z[0])); 
  q = p->back->next->next->back; 
  outerBranches.push_back(Branch(q->number, q->back->number, q->z[0])); 

  bool foundA = false; 
  bool foundB = false; 
  for(auto &v : outerBranches)
    {
      foundA |= v.nodeIsInBranch(switching.first); 
      foundB |= v.nodeIsInBranch(switching.second); 
    }
} 

void NniMove::multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<std::shared_ptr<AbstractPrior> > brPr) const 
{
  multiplyBranch(innerBranch, traln, hastings, prior, rand, multiplier,brPr, "someNNI"); 
  for(auto &v : outerBranches)
    multiplyBranch(v,traln, hastings, prior, rand, multiplier, brPr, "someNNI"); 
} 


void NniMove::multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter , std::vector<std::shared_ptr<AbstractPrior> > priors, std::string name) const 
{
  nodeptr p = branch.findNodePtr(traln); 

  double multi = rand.drawMultiplier(parameter); 
  double oldV = traln.getBranch(p).getLength(); 
  double newV = pow(oldV, multi); 

  traln.setBranchLengthBounded(newV, 0, p); 
  
  double realM = log(newV) /  log(oldV);   
  AbstractProposal::updateHastings(hastings, realM, name); 
  
  assert(priors.size() == 1); 
  auto brPr = priors[0]; 
  prior.updateBranchLengthPrior(traln,oldV, newV, brPr);
}


