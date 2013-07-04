#pragma once 
#include <vector>
#include "TreeAln.hpp"

// TODO this class is almost useless ... 


class Topology
{
public: 
  explicit Topology(int numTax); 
  Topology& operator=(const Topology& rhs); 
  
  void saveTopology(TreeAln &traln); 
  void restoreTopology(TreeAln &traln);   

private: 
  std::vector<Branch> branches;   
}; 
 
