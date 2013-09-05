#ifndef BIPARTITION_EXTRACTOR_HPP
#define BIPARTITION_EXTRACTOR_HPP

#include <unordered_map>

#include "TreeProcessor.hpp"
#include "BipartitionHash.hpp"

class BipartitionExtractor : public TreeProcessor
{
public: 
  BipartitionExtractor(std::vector<std::string> files);
  BipartitionExtractor( BipartitionExtractor&& rhs) ; 
  BipartitionExtractor& operator=(BipartitionExtractor rhs); 

  void extractBipsNew(); 
  void printBipartitions(std::string id) const ;
  void printBipartitionStatistics(std::string id) const ; 
  void printFileNames(std::string id) const ; 
  void printBranchLengths(std::string id) const ; 

  const std::vector<BipartitionHashNew>& getBipartitionHashes() const {return bipHashes; } 

  std::string bipartitionsToTreeString(std::vector<Bipartition> bips, bool printSupport) const ; 
  void buildTreeRecursive(nat currentId, const std::vector<std::vector<nat> >& directSubBips, const std::vector<Bipartition> &bips, std::stringstream &result, bool printSupport) const; 

private: 
  void extractUniqueBipartitions() ;
  nat getNumTreesInFile(std::string file) const ; 

private: 			// ATTRIBUTES
  std::vector<BipartitionHashNew> bipHashes; 
  std::unordered_map<Bipartition,nat> uniqueBips; 
}; 


#endif
