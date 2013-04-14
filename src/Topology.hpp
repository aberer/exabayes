#pragma once 
#include <vector>
#include "TreeAln.hpp"



/** 
    @brief stores a topology and the associated branch lengths 

 */
class Topology
{
public: 
  explicit Topology(int numTax); 
  ~Topology();
  
  void saveTopology(TreeAln &traln); 
  void restoreTopology(TreeAln &traln);   

private: 
  vector<branch*> branches;   

  void traverseAndSave(TreeAln &traln, nodeptr p, nat &number); 

}; 
 
