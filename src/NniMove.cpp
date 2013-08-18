#include "NniMove.hpp"
#include "LikelihoodEvaluator.hpp"
#include "TreePrinter.hpp"
#include "AbstractProposal.hpp"

void NniMove::applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const 
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

  traln.clipNode(pA,pBBack);  
  traln.setBranch(traln.getBranch(pA, params),params); 
  traln.clipNode(pB, pABack); 
  traln.setBranch(traln.getBranch(pB, params),params); 

} 

void NniMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const 
{
  applyToTree(traln, params);

  // revert branch lengths 
  // nodeptr p = innerBranch.findNodePtr(traln); 

  // reset inner branch 
  // traln.clipNode(p,p->back); 
  traln.setBranch(innerBranch, params); 

  for(auto &elem : outerBranches)
    {
      nodeptr q = elem.findNodePtr(traln); 
      traln.clipNode(q,q->back); 
      traln.setBranch(elem,params);
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



void NniMove::extractMoveInfo(const TreeAln &traln, std::vector<Branch> description, const std::vector<AbstractParameter*> &params) 
{
  assert(description.size() == 2 ) ; 
  std::pair<int,int>  _switching = { description.at(0).findNodePtr(traln)->next->back->number , description.at(1).findNodePtr(traln)->number}; 
  Branch _innerBranch = description.at(0); 

  outerBranches.clear(); 

  innerBranch = _innerBranch; 
  switching = _switching; 

  innerBranch = traln.getBranch(innerBranch.findNodePtr(traln), params);

  auto p = innerBranch.findNodePtr(traln); 

  for(auto &ptr : {p->next->back, p->next->next->back, 
	p->back->next->back, p->back->next->next->back })    
    {
      Branch b = traln.getBranch(Branch(ptr->number, ptr->back->number), params); 
      outerBranches.push_back(b);   
    }

  bool foundA = false; 
  bool foundB = false; 
  for(auto &v : outerBranches)
    {
      foundA |= v.nodeIsInBranch(switching.first); 
      foundB |= v.nodeIsInBranch(switching.second); 
    }
} 

#if 0 
void NniMove::multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, const std::vector<AbstractParameter*> &params) const 
{
  multiplyBranch(innerBranch, traln, hastings, prior, rand, multiplier,brPr, "someNNI"); 
  for(auto &v : outerBranches)
    multiplyBranch(v,traln, hastings, prior, rand, multiplier, brPr, "someNNI"); 
} 
#endif

#if 0 
void NniMove::multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, 
			     Randomness &rand, double parameter , std::vector<AbstractPrior* > priors, std::string name) const 
{
  nodeptr p = branch.findNodePtr(traln); 

  double multi = rand.drawMultiplier(parameter); 
  double oldV = traln.getBranch(p).getLength(); 
  double newV = pow(oldV, multi); 

  traln.setBranch(Branch(p->number, p->back->number, newV));
  
  double realM = log(newV) /  log(oldV);   
  AbstractProposal::updateHastings(hastings, realM, name); 
  
  assert(priors.size() == 1); 
  auto brPr = priors[0]; 
  prior.updateBranchLengthPrior(traln,oldV, newV, brPr);
}
#endif




NniMove NniMove::getInvertseMove() const 
{
  return NniMove(*this); 
} 
