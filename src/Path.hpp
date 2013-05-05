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

  bool nodeIsOnPath(int node);  
  void saveBranchLengthsPath(TreeAln& traln); 
  bool isOuterNode(int node); 
  void debug_assertPathExists(TreeAln& traln); 
  void pushToStackIfNovel(branch b, int numTip); 
  void clearStack(); 
  void pushStack(branch value); 
  branch popStack(); 
  bool stackIsEmpty(); 
  int size() {return stack.size(); }
  branch& peekStack(); 
  void restoreBranchLengthsPath(TreeAln &traln); 

  int getNthNodeInPath(nat num) ; 
  int getNumberOfNodes() const {return stack.size()  + 1 ;   }

  void printWithBLs(TreeAln &traln ); 
  
  // TDOO avoid that later 
  // branch& operator[](int num){return stack[num];  }
  branch& at(int num){return stack[num]; }

  // TODO this should return a new path 
  void multiplyBranch(TreeAln &traln, Randomness &rand, branch b, double parameter, double &hastings); 

  friend ostream& operator<<(ostream &out, const Path &rhs)  ;
  void destroyOrientationAlongPath(tree *tr,  nodeptr p); 

  // TODO much of that should not be here -- most of it is, because of conversion 

  // ugly 
  vector<branch>& getStack(){return stack; }

private: 
  vector<branch> stack; 

}; 


#endif
