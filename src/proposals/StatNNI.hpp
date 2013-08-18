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

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) ; 
  virtual void resetState(TreeAln &traln) ; 

  // virtual Branch prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return Branch(0,0);}
  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}

  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  
  
  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out) const { } //disabled


private: 			// METHODS
  void treatOneBranch(nodeptr p, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand); 

private:			// ATTRIBUTES
  double multiplier; 
  Path path; 
  NniMove move; 
#ifdef _EXPERIMENTAL_INTEGRATION_MODE
  nat invocations; 
#endif
};

#endif
