#ifndef _TREE_PRINTER_HPP
#define _TREE_PRINTER_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include <sstream>


class TreePrinter
{
public: 
  TreePrinter(bool withBranchLengths , bool withInternalNodes , bool withRealNames) 
    : withBranchLengths(withBranchLengths), withInternalNodes(withInternalNodes), withRealNames(withRealNames){}

  string printTree(const TreeAln &traln ); 

private: 
  bool withBranchLengths; 
  bool withInternalNodes; 
  bool withRealNames; 

  void helper(const TreeAln &traln, stringstream &ss, nodeptr p); 
  

}; 





#endif
