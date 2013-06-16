#ifndef _NNI_MOVE
#define _NNI_MOVE

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "PriorBelief.hpp"

class  NniMove
{
public: 
  void apply(TreeAln &traln) const ; 
  void revert(TreeAln &traln) const; 
  void disortient(TreeAln &traln) const ; 
  void init(const TreeAln &traln, const Branch& _innerBranch, pair<int,int> _switching); 
  Branch getEvalBranch() const {return innerBranch; }
  
  void multiplyAllBranches( TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter, vector<shared_ptr<AbstractPrior> > priors, string name) ; 

  
private: 
  Branch innerBranch; 
  vector<Branch> outerBranches; 
  pair<int,int> switching; 


  void multiplyBranch(const Branch &branch, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand, double parameter , vector<shared_ptr<AbstractPrior> > priors, string name); 
}; 



#endif
