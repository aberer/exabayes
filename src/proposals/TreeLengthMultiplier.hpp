#include "AbstractProposal.hpp"


class TreeLengthMultiplier : public AbstractProposal
{
public: 
  TreeLengthMultiplier( double _multiplier) ; 
  virtual ~TreeLengthMultiplier(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior)  ; 
  virtual void autotune() ;

  virtual AbstractProposal* clone() const;  

  virtual void readFromCheckpointCore(std::istream &in);
  virtual void writeToCheckpointCore(std::ostream &out) const;

  // virtual Branch prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return Branch(0,0);}
  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}

private: 			// METHODS
  void multiplyBranchLengthsRecursively(TreeAln& traln, nodeptr p, double multiHere); 

private: 			// ATTRIBUTES
  double multiplier; 		// the tuning variable  
  double initTreeLength; 
  std::vector<Branch> storedBranches; 
} ; 
