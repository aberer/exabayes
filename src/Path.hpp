/** 
    @file Path.hpp represents a  path in a tree
    
    Currently we are assuming, the path always is connected. 
 */ 


#ifndef _PATH_H
#define _PATH_H

#include <vector>

#include "axml.h"
#include "branch.h"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "PriorBelief.hpp"

class Path
{
public: 
  // lets do this correctly for once ... 
  Path();
  ~Path();   
  Path& operator=(const Path &rhs ) ;   
  Path( const Path &rhs) ; 
  friend void swap(Path &first, Path &second); 
  
/** @brief returns true, if the node with a given id is part of this branch */ 
  bool nodeIsOnPath(int node) const;  
  /** @brief for all branches in the path, copy over the branch lengths */ 
  void saveBranchLengthsPath(const TreeAln& traln); 

  /** @brief asserts that this path exists in a given tree */ 
  void debug_assertPathExists(TreeAln& traln); 

  /** @brief assigns stored branch lengths of a path to a given tree  */ 
  void restoreBranchLengthsPath(TreeAln &traln ,PriorBelief &prior) const ; 

  /** @brief only add a branch to the path, if it is novel. If the new
      branch cancels out an existing branch, the path is shortened again */ 
  void pushToStackIfNovel(branch b, int numTip); 

  // straight-forward container methods 
  void append(branch value); 
  void clear(); 
  /** @brief number of branches in the path */ 
  nat size() const {return stack.size(); }

  /** @brief yields the branch */  
  branch& at(int num){return stack[num]; }
  branch at(int num) const{return stack[num];}

  /** @brief reverse the path */ 
  void reverse(); 

  /** @brief removes the last element */ 
  void pop(); 
  /** @brief removes the first element */ 
  void popFront(); 
  

  /** @brief returns the id of the nth node in the path. nodes 0 and n+1 are the outer nodes in this path that do not have a neighbor. */ 
  int getNthNodeInPath(nat num) const ; 
  
  /** @brief gets the number of nodes represented by the path (assuming it is connected)  */
  int getNumberOfNodes() const {return stack.size()  + 1 ;   }
  void printWithBLs(TreeAln &traln ) const; 

  void multiplyBranch(TreeAln &traln, Randomness &rand, branch b, double parameter, double &hastings, PriorBelief &prior, shared_ptr<AbstractPrior> prBr) const; 

  // // TODO this should return a new path instance  
  // void multiplyBranch(TreeAln &traln, Randomness &rand, branch b, double parameter, double &hastings, PriorBelief &prior); 

  void findPath(const TreeAln& traln, nodeptr p, nodeptr q);

  friend ostream& operator<<(ostream &out, const Path &rhs)  ;

private: 
  vector<branch> stack; 

  bool findPathHelper(const TreeAln &traln, nodeptr p, const branch &target);

}; 


#endif


