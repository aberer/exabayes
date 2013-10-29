#ifndef _TREE_PROCESSOR_HPP
#define _TREE_PROCESSOR_HPP

#include <string>
#include "TreeAln.hpp"

class TreeProcessor
{
public: 
  TreeProcessor(std::vector<std::string> fileNames); 
  TreeProcessor(TreeProcessor&& tp) ; 
  TreeProcessor& operator=(TreeProcessor &&tp); 
  
  const std::vector<std::string> getTaxa() const {return taxa; }

protected: 			// METHODS
  void fillTaxaInfo(std::string fileName); 
  template<bool readBl>
  void nextTree(std::istream &treefile); 
  void skipTree(std::istream &iss); 
  /** 
      @brief only initializes a raw tree, no partitions or alignment information. 
      important: does NOT need a bytefile 
  */ 
  void initializeTreeOnly(int numTax);  
  static std::string trim(const std::string& str, const std::string& whitespace  = " \t"); 

protected: 			// ATTRIBUTES
  std::unique_ptr<TreeAln> tralnPtr;
  std::vector<std::string> fns; 
  std::vector<std::string> taxa; 
}; 


#endif
