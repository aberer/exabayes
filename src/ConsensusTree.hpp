#ifndef CONSENSUSTREE_HPP
#define CONSENSUSTREE_HPP

#include "BipartitionExtractor.hpp"
#include <vector>
#include <string>

class ConsensusTree 
{
public: 
  ConsensusTree(std::string file); 
  std::string getConsensusTreeString(double threshold, bool isMRE); 
  std::vector<Bipartition>  getRefinedConsensusTree(const std::vector<Bipartition> &consensusBips, const std::vector<Bipartition> &minorityBips) const ; 

private:   
  BipartitionExtractor bipEx; 
  nat totalTrees; 
}; 


#endif
