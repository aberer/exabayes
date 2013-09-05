#include "NniMove.hpp"
#include "LikelihoodEvaluator.hpp"
#include "TreePrinter.hpp"
#include "AbstractProposal.hpp"

void NniMove::applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const 
{
  auto a = BranchPlain(switching.first, innerBranch.getPrimNode()); 
  auto b = BranchPlain{}; 
  if( not a.exists(traln))
    {
      a = BranchPlain(switching.first, innerBranch.getSecNode()); 
      assert(a.exists(traln)); 
      b = BranchPlain(switching.second, innerBranch.getPrimNode()); 
      assert(b.exists(traln));       
    }
  else 
    {
      b = BranchPlain(switching.second, innerBranch.getSecNode()); 
      assert(b.exists(traln));       
    }  

  nodeptr pA = a.findNodePtr(traln),
    pB = b.findNodePtr(traln), 
    pABack = pA->back,
    pBBack = pB->back; 
  
  traln.clipNode(pA,pBBack);  
  for(auto &param : params)
    {
      auto elem = traln.getBranch(pA, param); 
      traln.setBranch(elem, param); 
    }
  
  traln.clipNode(pB, pABack); 
  for(auto &param : params)
    {
      auto elem = traln.getBranch(pB, param); 
      traln.setBranch(elem, param); 
    }
} 

void NniMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const 
{
  applyToTree(traln, params);

  // reset inner branch 
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
  if( BranchPlain(p->number, p->back->number).equalsUndirected(innerBranch.toPlain()))
    {      
      disorientHelper(traln, p); 
      disorientHelper(traln, p->back); 
    }
  else if(  innerBranch.hasNode(p->number) )
    {
      assert(not innerBranch.hasNode( p->back->number)); 
      disorientAtNode(traln, p);
      int node = innerBranch.getOtherNode(p->number); 
      disorientAtNode(traln, BranchPlain(node, p->number).findNodePtr(traln)); 
    } 
  else 
    {
      // disallowing evaluation on any more peripheral pointers. I guess this would not be needed.   
      assert(0); 
    }
} 



std::ostream& operator<<(std::ostream& out, const NniMove &move)
{
  // weird errors 
#if 0 
  // out << move.innerBranch.toString(); 
  auto bl = move.innerBranch.toPlain(); 
  out << bl; 

  out << ";" ; 
  bool isFirst = true; 
  for(auto & b : move.outerBranches)
    {
      if(isFirst)
	isFirst = false; 
      else 
	out << ","; 
      out << b.toPlain();
    }
  
  out << ";"; 
  out << move.switching.first << "," << move.switching.second; 
#endif 
  return out; 
}



void NniMove::extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description, const std::vector<AbstractParameter*> &params) 
{
  assert(description.size() == 2 ) ; 
  std::pair<int,int>  _switching = { description.at(0).findNodePtr(traln)->next->back->number , description.at(1).findNodePtr(traln)->number}; 
  auto _innerBranch = description.at(0);
					     
  outerBranches.clear(); 

  assert(params.size() == 1 ); 

  switching = _switching; 
  innerBranch = traln.getBranch(_innerBranch, params);

  auto p = innerBranch.findNodePtr(traln); 

  for(auto &ptr : {p->next->back, p->next->next->back, 
	p->back->next->back, p->back->next->next->back })    
    {
      auto b = traln.getBranch(BranchPlain(ptr->number, ptr->back->number), params); 
      outerBranches.push_back(b);   
    }

  bool foundA = false; 
  bool foundB = false; 
  for(auto &v : outerBranches)
    {
      foundA |= v.hasNode(switching.first); 
      foundB |= v.hasNode(switching.second); 
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



template<typename BT>
BranchPlain NniMove::mapToBranchBeforeMove(const BT &b) const 
{
  if( ( b.hasNode(switching.first) || b.hasNode(switching.second)  ) 
      && ( b.hasNode(innerBranch.getPrimNode()) || b.hasNode(innerBranch.getSecNode())) )
    {
      auto firstNode = switching.first; 
      if(not b.hasNode(firstNode))
	firstNode = switching.second; 
      auto otherNode = b.getOtherNode(firstNode) ; 
      auto newOtherNode = innerBranch.getOtherNode(otherNode); 
      return BranchPlain(firstNode, newOtherNode); 
    }
  else 
    return b; 
}



