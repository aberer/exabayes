#include <iostream>
#include <fstream>

#include "Path.hpp"
#include "mcmc/Chain.hpp"
#include "proposals/AbstractProposal.hpp"


void Path::clear()
{
  stack.clear(); 
}

void Path::append(BranchPlain value)
{
  stack.push_back(value);
}



void Path::pop()
{
  stack.pop_back();
}

void Path::popFront()
{
  stack.erase(stack.begin()); 
}


void Path::pushToStackIfNovel(BranchPlain b, const TreeAln &traln)
{  
  assert(stack.size() >= 2 ); 

  auto bPrev = at(size()-1); 
  if(stack.size() == 2)
    {
      if(not b.equalsUndirected(bPrev))
	append(b); 
    }
  else 
    {      
      if(b.equalsUndirected(bPrev))
	stack.pop_back();
      else if(bPrev.isTipBranch(traln))
	{
	  stack.pop_back();
	  append(b); 
	}
      else 
	append(b ); 
    }
}




void Path::debug_assertPathExists(TreeAln& traln)
{
#ifdef DEBUG_CHECK_TREE_CONSISTENCY
  tree *tr = traln.getTr();
  for(auto b : stack)
    assert(branchExists(tr, b)); 
#endif
}




/**
    @brief saves all branch lengths along the path in s. 
 */ 
void Path::saveBranchLengthsPath(const TreeAln& traln, const std::vector<AbstractParameter*> &params)
{
  bls.clear();
  for(auto &b : stack)
    {
      auto p = b.findNodePtr(traln);
      auto bl = traln.getBranch(p,params); 
      // tout << "saving " << bl << std::endl; 
      bls.push_back(bl); 
    }
}


bool Path::nodeIsOnPath(int node) const
{
  for(auto b : stack)
    if(b.hasNode(node))
      return true;       

  return false; 
}


std::ostream& operator<<(std::ostream &out, const Path &rhs)  
{
  for(auto b : rhs.stack)
    out << "(" << b.getPrimNode() << "," << b.getSecNode() << "),";       
  return out; 
}


void Path::multiplyBranch(TreeAln &traln, Randomness &rand, BranchLength b, double parameter, double &hastings, PriorBelief &prior, AbstractParameter* const param) const 
{  
#if 0 
  nodeptr p = b.findNodePtr(traln); 
  double multiplier = rand.drawMultiplier(parameter); 

  double oldZ = traln.getBranch(p, param).getLength(param); 
  
  double newZ = multiplier * oldZ; 

  assert(0);
  // , newZ
  traln.clipNode(p,p->back);   

  prior.updateBranchLengthPrior(traln, oldZ, newZ, param);

  double realMultiplier = log(newZ) / log(oldZ);    
  AbstractProposal::updateHastings(hastings, realMultiplier, "pathMod");; 
#else 
  assert(0); 
#endif
}


void Path::restoreBranchLengthsPath(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const 
{
  for(auto &bl : bls )
    {
      // tout << "restore " << bl << std::endl; 
      traln.setBranch(bl,blParams); 
    }
}


int Path::getNthNodeInPath(nat num)   const
{ 
  int result; 
  
  assert(num <= stack.size( )+  1 ); 
  // TODO stronger warrenty when constructing 
  // assert(stack.size() != 1 ); 
  
  if(num == 0)
    {
      auto  b = stack[0]; 
      if(stack[1].hasNode(b.getPrimNode()))
	result= b.getSecNode(); 
      else 
	{
	  assert(stack[1].hasNode(b.getSecNode())); 
	  result= b.getPrimNode() ; 
	}
    }
  else if(num == stack.size() )
    {
      auto b = stack[num-1]; 
      
      if(stack[num-2].hasNode(b.getPrimNode() ))
	result=  b.getSecNode(); 
      else 
	{
	  assert(stack[num-2].hasNode(b.getSecNode() )); 
	  result= b.getPrimNode(); 
	}
    }
  else 
    {      
      auto b = stack[num]; 
      if(stack[num-1].hasNode(b.getPrimNode() ))
	result= b.getPrimNode(); 
      else 
	{
	  assert(stack[num-1].hasNode(b.getSecNode() )); 
	  result= b.getSecNode(); 
	}	
    }

  return result; 
}


void Path::reverse()
{
  std::reverse(stack.begin(), stack.end());
}




bool Path::findPathHelper(const TreeAln &traln, nodeptr p, const BranchPlain &target)
{
  auto  curBranch = BranchPlain(p->number, p->back->number); 
  if( curBranch.equalsUndirected( target) ) 
    return true; 
  
  bool found = false; 
  for(auto q = p->next ; p != q && not found ; q = q->next)
    {
      found = findPathHelper(traln, q->back, target); 
      if(found)	
	{
	  append(BranchPlain(q->number, q->back->number));
	  return found; 
	}
    }
  return false; 
}



void Path::findPath(const TreeAln& traln, nodeptr p, nodeptr q)
{  
  // cout << "trying to find path between " << p->number << "/" << p->back->number << " and "  << q->number << "/" << q->back->number  << endl; 
  
  stack.clear();   

  auto targetBranch = BranchPlain(q->number, q->back->number ); 

  bool found = findPathHelper(traln, p->back, targetBranch); 
  if(found)  
    {      
      append(BranchPlain(p->number, p->back->number)); 
      return ; 
    }

  assert(stack.size() == 0); 
  
  found = findPathHelper(traln, p->next->back, targetBranch); 
  if(found  )
    {
      append(BranchPlain(p->next->number, p->next->back->number)); 
      return; 
    }

  assert(stack.size() == 0);   
  found = findPathHelper(traln, p->next->next->back, targetBranch); 
  if(found)
    append(BranchPlain(p->next->next->number, p->next->next->back->number)); 
  assert(found);
}
