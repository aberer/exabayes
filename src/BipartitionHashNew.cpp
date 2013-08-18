#include "BipartitionHashNew.hpp"
#include <algorithm>
#include <iostream>

#include "Branch.hpp"
#include "Randomness.hpp"

BipartitionHashNew::BipartitionHashNew(nat numTax)   
  :bipMeaning(numTax)
  , treesAdded(0)
{
  // must be deterministic!
  auto seed = randCtr_t{0xdeadbeef,0xbadcafee}; 
  Randomness rand(seed);

  // we do nat allow a zero 
  for(auto &r : bipMeaning)
    while(r == 0)
      r = rand();
}


void BipartitionHashNew::addTree(const TreeAln &traln, bool withBranch)
{
  auto pStart = traln.getTr()->nodep[1]->back; 
  auto desc = traln.getDescendents(Branch(pStart->number, pStart->back->number)); 

  addElement(traln, desc.first.getInverted().findNodePtr(traln), withBranch); 
  addElement(traln, desc.second.getInverted().findNodePtr(traln), withBranch); 

  ++treesAdded;
}


Bipartition BipartitionHashNew::addElement(const TreeAln &traln, nodeptr p, bool withBranch)
{
  auto curBranch = Branch(p->number, p->back->number);

  if(curBranch.isTipBranch(traln))
    {
      auto result = Bipartition{}; 
      // TODO maybe insert trivial bipartitions as well 
      result.initializeWithTaxon(p->number-1, bipMeaning[p->number-1]); 
      return result;  
    }
  
  auto desc = traln.getDescendents(curBranch);

  auto bipA = addElement(traln, desc.first.getInverted().findNodePtr(traln), withBranch); 
  auto bipB = addElement(traln, desc.second.getInverted().findNodePtr(traln), withBranch) ; 

  auto result = bipA | bipB; 

  if(withBranch)
    {
      // assert(0); 
      assert(traln.getNumBranches() == 1); 
      // notice: if we read, we only want to have one branch length
      // there => non-issue, but be careful
      auto bl = p->z[0]; 	// TODO 
      bipBranchLengths[result].push_back(bl); 
    }
  
  auto &precBip = bipPresence[result]; 
  precBip.reserve(treesAdded);
  precBip.set(treesAdded);

  return result; 
} 


std::vector<double>
BipartitionHashNew::getBranchLengths(const Bipartition& bip) const
{
  auto result = std::vector<double>{}; 
  auto iter = bipBranchLengths.find(bip); 
  if(iter != bipBranchLengths.end())
    result = iter->second; 
  return result;
}
