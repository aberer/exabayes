#include <iostream>
#include <fstream>

#include "Path.hpp"
#include "Chain.hpp"
#include "AbstractProposal.hpp"




void Path::clear()
{
  stack.clear(); 
}

void Path::append(Branch value)
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


void Path::pushToStackIfNovel(Branch b, const TreeAln &traln)
{  
  assert(stack.size() >= 2 ); 

  Branch bPrev = at(size()-1); 
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
void Path::saveBranchLengthsPath(const TreeAln& traln)
{
  auto *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches() ;
  assert(numBranches == 1 ); 

  for(auto& b : stack)
    {
      nodeptr p = b.findNodePtr(traln); 

      // for(int j = 0; j < numBranches; ++j)
      // 	{
	  double tmp = traln.getBranchLength( p,0); 
	  b.setLength( tmp) ; 
	// }
    }
}


bool Path::nodeIsOnPath(int node) const
{
  for(auto b : stack)
    if(b.nodeIsInBranch(node))
      return true;       

  return false; 
}


ostream& operator<<(ostream &out, const Path &rhs)  
{
  for(auto b : rhs.stack)
    out << "(" << b.getPrimNode() << "," << b.getSecNode() << "),";       
  return out; 
}



void Path::multiplyBranch(TreeAln &traln, Randomness &rand, Branch b, double parameter, double &hastings, PriorBelief &prior, PriorPtr brPr) const 
{  
  tree *tr = traln.getTr(); 
  nodeptr p = b.findNodePtr(traln); 
  double multiplier = rand.drawMultiplier(parameter); 

  double oldZ = traln.getBranchLength(p,0);
  double newZ = multiplier * oldZ; 

  traln.clipNode(p,p->back, newZ);   

  prior.updateBranchLengthPrior(traln, oldZ, newZ, brPr);

  double realMultiplier = log(newZ) / log(oldZ);    
  AbstractProposal::updateHastings(hastings, realMultiplier, "pathMod");; 
}


void Path::restoreBranchLengthsPath(TreeAln &traln ,PriorBelief &prior) const 
{
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 
  for(auto b : stack)
    {
      nodeptr p = b.findNodePtr(traln); 
      double tmp = b.getLength(); 
      traln.clipNode(p, p->back, tmp); 
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
      Branch b = stack[0]; 
      if(stack[1].nodeIsInBranch(b.getPrimNode()))
	result= b.getSecNode(); 
      else 
	{
	  assert(stack[1].nodeIsInBranch(b.getSecNode())); 
	  result= b.getPrimNode() ; 
	}
    }
  else if(num == stack.size() )
    {
      Branch b = stack[num-1]; 
      
      if(stack[num-2].nodeIsInBranch(b.getPrimNode() ))
	result=  b.getSecNode(); 
      else 
	{
	  assert(stack[num-2].nodeIsInBranch(b.getSecNode() )); 
	  result= b.getPrimNode(); 
	}
    }
  else 
    {      
      Branch b = stack[num]; 
      if(stack[num-1].nodeIsInBranch(b.getPrimNode() ))
	result= b.getPrimNode(); 
      else 
	{
	  assert(stack[num-1].nodeIsInBranch(b.getSecNode() )); 
	  result= b.getSecNode(); 
	}	
    }

  return result; 
}


void Path::reverse()
{
  std::reverse(stack.begin(), stack.end());
}




bool Path::findPathHelper(const TreeAln &traln, nodeptr p, const Branch &target)
{
  Branch curBranch =  Branch(p->number, p->back->number); 
  if( curBranch.equalsUndirected( target) ) 
    return true; 
  
  bool found = false; 
  for(auto q = p->next ; p != q && not found ; q = q->next)
    {
      found = findPathHelper(traln, q->back, target); 
      if(found)	
	{
	  append(Branch(q->number, q->back->number));
	  return found; 
	}
    }
  return false; 
}



void Path::findPath(const TreeAln& traln, nodeptr p, nodeptr q)
{  
  // cout << "trying to find path between " << p->number << "/" << p->back->number << " and "  << q->number << "/" << q->back->number  << endl; 
  
  stack.clear();   

  Branch targetBranch = Branch(q->number, q->back->number ); 

  bool found = findPathHelper(traln, p->back, targetBranch); 
  if(found)  
    {      
      append(Branch(p->number, p->back->number)); 
      return ; 
    }

  assert(stack.size() == 0); 
  
  found = findPathHelper(traln, p->next->back, targetBranch); 
  if(found  )
    {
      append(Branch(p->next->number, p->next->back->number)); 
      return; 
    }

  assert(stack.size() == 0);   
  found = findPathHelper(traln, p->next->next->back, targetBranch); 
  if(found)
    append(Branch(p->next->next->number, p->next->next->back->number)); 
  assert(found);
}
