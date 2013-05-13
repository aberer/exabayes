#include <iostream>
#include <fstream>

#include "Path.hpp"
#include "branch.h"
#include "Chain.hpp"


Path::Path()
{
}

Path::~Path()
{
}


Path::Path(const Path &rhs)
    : stack(rhs.stack)
{    
}



Path& Path::operator=(const Path &rhs )  
{
  Path tmp(rhs); 
  swap(*this, tmp); 
  return *this; 
}





void Path::clear()
{
  stack.clear(); 
}

void Path::append(branch value)
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


void Path::pushToStackIfNovel(branch b, int numTip)
{  
  assert(stack.size() >= 2 ); 

  branch bPrev = at(size()-1); 
  if(stack.size() == 2)
    {
      if(NOT branchEqualUndirected(b, bPrev))
	append(b); 
    }
  else 
    {      
      if(branchEqualUndirected(b,bPrev))
	stack.pop_back();
      else if(isTipBranch(bPrev, numTip))
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
void Path::saveBranchLengthsPath(TreeAln& traln)
{
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches() ;

  for(branch& b : stack)
    {
      nodeptr p = findNodeFromBranch(tr, b); 

      for(int j = 0; j < numBranches; ++j)
	{
	  double tmp = traln.getBranchLength( p,0); 
	  b.length[j] = tmp ; 
	}
    }
}


bool Path::nodeIsOnPath(int node)
{
  for(auto b : stack)
    if(nodeIsInBranch(node, b))
      return true;       

  return false; 
}


ostream& operator<<(ostream &out, const Path &rhs)  
{
  for(auto b : rhs.stack)
    out << "(" << b.thisNode << "," << b.thatNode << "),";       
  return out; 
}



void Path::multiplyBranch(TreeAln &traln, Randomness &rand, branch b, double parameter, double &hastings, PriorBelief &prior)
{  
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches();
  nodeptr  p = findNodeFromBranch(tr, b); 
  double multiplier = rand.drawMultiplier(parameter); 

  double oldZ = branchLengthToReal(tr, traln.getBranchLength( p,0)); 
  double newZ = branchLengthToInternal(tr, multiplier * oldZ); 
#ifdef PRINT_MULT
  cout  << setprecision(6) << "spr: " << oldZ <<   " * " << multiplier << " = "  << multiplier * oldZ << endl; // 
#endif

  assert(not prior.believingInFixedBranchLengths()); 

  prior.updateBranchLength(oldZ, branchLengthToReal(traln.getTr(),newZ));
  
  updateHastings(hastings, multiplier, "pathMod");; 
  hookup(p,p->back, &newZ, numBranches);   
}


void Path::restoreBranchLengthsPath(TreeAln &traln ,PriorBelief &prior)
{
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 
  for(auto b : stack)
    {
      nodeptr p = findNodeFromBranch(traln.getTr(), b); 
      if(not prior.believingInFixedBranchLengths())
	prior.updateBranchLength(branchLengthToReal(traln.getTr(), p->z[0]), branchLengthToReal(traln.getTr(), b.length[0]) );
      hookup(p, p->back, b.length, numBranches); 
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
      branch b = stack[0]; 
      if(nodeIsInBranch(b.thisNode, stack[1]))
	result= b.thatNode; 
      else 
	{
	  nodeIsInBranch(b.thatNode,stack[1]); 
	  result= b.thisNode; 
	}
    }
  else if(num == stack.size() )
    {
      branch b = stack[num-1]; 
      
      if(nodeIsInBranch(b.thisNode, stack[num-2]))
	result=  b.thatNode; 
      else 
	{
	  assert(nodeIsInBranch(b.thatNode, stack[num-2])); 
	  result= b.thisNode; 
	}
    }
  else 
    {
      branch b = stack[num]; 
      if(nodeIsInBranch(b.thisNode, stack[num-1]))
	result= b.thisNode; 
      else 
	{
	  assert(nodeIsInBranch(b.thatNode, stack[num-1])); 
	  result= b.thatNode; 
	}	
    }

  return result; 
}

void swap(Path &first, Path &second)
{
  using std::swap; 
  swap(first.stack, second.stack); 
}


void Path::reverse()
{
  std::reverse(stack.begin(), stack.end());
}
