#ifndef BIPARTITIONHASHNEW_H
#define BIPARTITIONHASHNEW_H

#include <unordered_map>

#include "axml.h"
#include "TreeAln.hpp"
#include "Bipartition.hpp"


class BipartitionHash
{
public: 
  BipartitionHash(nat numTax); 
  void addTree(const TreeAln &traln, bool withBranch, bool withTrivial); 

  std::unordered_map<Bipartition, Bipartition>::const_iterator begin() const{return bipPresence.begin(); }
  std::unordered_map<Bipartition, Bipartition>::const_iterator end() const{return bipPresence.end(); }

  std::vector<double> getBranchLengths(const Bipartition& bip) const;
  Bipartition getPresence(const Bipartition &bip) const; 

  nat getTreesAdded() const {return treesAdded; }

private: 			// METHODS
  Bipartition addElement(const TreeAln &traln, nodeptr p, bool withBranch, bool withTrivial); 

private: 			// ATTRIBUTES
  std::unordered_map<Bipartition, Bipartition> bipPresence; 
  std::unordered_map<Bipartition, std::vector<double> > bipBranchLengths; 
  std::vector<nat> bipMeaning; 
  nat treesAdded; 
}; 


#endif
