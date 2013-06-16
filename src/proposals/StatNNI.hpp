/** 
    @brief our flavour of the statistical nearest neighbor interchange

    in constrast to mb, this proposal will always induce a topological
    change
    
    @notice we could multiply some additional branches
 */ 


#ifndef _STAT_NNI_H
#define _STAT_NNI_H

class Chain; 
#include "TreeAln.hpp"
#include "Path.hpp"
#include "AbstractProposal.hpp"
#include "NniMove.hpp"

class StatNNI : public AbstractProposal
{
public: 
  StatNNI( double multiplier);
  virtual ~StatNNI(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  

private:
  double multiplier; 
  Path path; 
  NniMove move; 

  void treatOneBranch(nodeptr p, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand); 
};

#endif