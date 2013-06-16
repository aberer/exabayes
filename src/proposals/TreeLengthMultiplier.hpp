#include "AbstractProposal.hpp"


class TreeLengthMultiplier : public AbstractProposal
{
public: 
  TreeLengthMultiplier( double _multiplier) ; 
  virtual ~TreeLengthMultiplier(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior)  ; 
  virtual void autotune() ;

  virtual AbstractProposal* clone() const;  

private: 
  double multiplier; 		// the tuning variable  
  double rememMultiplier; 	// for resetting 
  double initTreeLength; 
  
  void multiplyBranchLengthsRecursively(TreeAln& traln, nodeptr p, double multiHere); 

  vector<branch> storedBranches; 

} ; 
