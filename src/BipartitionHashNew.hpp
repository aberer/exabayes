#ifndef BIPARTITIONHASHNEW_H
#define BIPARTITIONHASHNEW_H

#include <unordered_map>

#include "axml.h"
#include "TreeAln.hpp"
#include "Bipartition.hpp"

// template by branch length 

class BipartitionHashNew
{
public: 
  BipartitionHashNew(nat numTax); 
  void addTree(const TreeAln &traln, bool withBranch); 

  std::unordered_map<Bipartition, Bipartition>::const_iterator begin() const{return bipPresence.begin(); }
  std::unordered_map<Bipartition, Bipartition>::const_iterator end() const{return bipPresence.end(); }

  std::vector<double> getBranchLengths(const Bipartition& bip) const;
  
private: 			// METHODS
  Bipartition addElement(const TreeAln &traln, nodeptr p, bool withBranch); 

private: 			// ATTRIBUTES
  std::unordered_map<Bipartition, Bipartition> bipPresence; 
  std::unordered_map<Bipartition, std::vector<double> > bipBranchLengths; 
  std::vector<nat> bipMeaning; 
  nat treesAdded; 
}; 


#endif
