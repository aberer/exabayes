/**
   @file branch.h
   
   @brief Everything that deals with branches.  
   

   TODO call by value for branches is too expensive  
 */


#ifndef _BRANCH_H
#define _BRANCH_H

#include <vector>
#include <iostream>
using namespace std; 

#include "axml.h"

class TreeAln; 

typedef struct  _branch
{
  int thisNode;  /// in case we want to give the branch some orientation, this is the node we want to extract   
  int thatNode;  
  
  double length[NUM_BRANCHES]; 	/* TODO expensive?  */
} branch; /// represents a branch 


void extractBranches(const TreeAln &traln, vector<branch> &result) ; 

branch invertBranch(branch b); 
branch constructBranch(int thisNode, int thatNode); 
boolean branchExists(tree *tr, branch b); 
nodeptr findNodeFromBranch(tree *tr, branch b ); 
branch findRoot(tree *tr); 
boolean isTipBranch(branch b, int numTip); 
boolean branchEqualUndirected(branch b1, branch b2); 
boolean nodeIsInBranch(int number, branch b ); 
int getIntersectingNode(branch b1, branch b2); 
branch getThirdBranch(tree *tr, branch b1, branch b2); 
int getOtherNode(int node, branch b); 
double branchLengthToInternal(tree *tr, double realBL);
double branchLengthToReal(tree *tr, double internalBL); 
void divideBranchLengthsWithRatio(tree *tr, double orig,  double ratio, double *resultA, double *resultB); 
double combineBranchLengths(tree *tr, double origA, double origB); 
double getRatio(tree *tr, double a, double b ); 

ostream& operator<<(ostream& rhs, const branch &b);

#endif
