/** 
    @file Path.hpp represents a  path in a tree
    
    Currently we are assuming, the path always is connected. 
 */ 


#ifndef _PATH_H
#define _PATH_H

#include <vector>

#include "axml.h"
#include "bayes.h"
#include "branch.h"
#include "TreeAln.hpp"

class Path
{
public: 
  Path();
  ~Path();   
  
  /** @brief returns true, if the node with a given id is part of this branch */ 
  bool nodeIsOnPath(int node);  
  /** @brief for all branches in the path, copy over the branch lengths */ 
  void saveBranchLengthsPath(TreeAln& traln); 

  /** @brief asserts that this path exists in a given tree */ 
  void debug_assertPathExists(TreeAln& traln); 

  /** @brief assigns stored branch lengths of a path to a given tree  */ 
  void restoreBranchLengthsPath(TreeAln &traln); 

  /** @brief only add a branch to the path, if it is novel. If the new
      branch cancels out an existing branch, the path is shortened again */ 
  void pushToStackIfNovel(branch b, int numTip); 

  // straight-forward container methods 
  void append(branch value); 
  void clear(); 
  int size() {return stack.size(); }
  branch& at(int num){return stack[num]; }

  /** @brief returns the id of the nth node in the path. nodes 0 and n+1 are the outer nodes in this path that do not have a neighbor. */ 
  int getNthNodeInPath(nat num) ; 
  
  /** @brief gets the number of nodes represented by the path (assuming it is connected)  */
  int getNumberOfNodes() const {return stack.size()  + 1 ;   }
  void printWithBLs(TreeAln &traln ); 

  // TODO this should return a new path instance  
  void multiplyBranch(TreeAln &traln, Randomness &rand, branch b, double parameter, double &hastings); 

  friend ostream& operator<<(ostream &out, const Path &rhs)  ;

  
  void destroyOrientationAlongPath(tree *tr,  nodeptr p); // TODO move? 


private: 
  vector<branch> stack; 

}; 


#endif
