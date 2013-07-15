#include "AbstractProposal.hpp"


class TreeLengthMultiplier : public AbstractProposal
{
public: 
  TreeLengthMultiplier( double _multiplier) ; 
  virtual ~TreeLengthMultiplier(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior)  ; 
  virtual void autotune() ;

  virtual AbstractProposal* clone() const;  

  virtual void readFromCheckpointCore(std::ifstream &in) { in >> multiplier;  } 
  virtual void writeToCheckpointCore(std::ofstream &out)const {out << multiplier << DELIM; } 

private: 
  double multiplier; 		// the tuning variable  
  double initTreeLength; 
  
  void multiplyBranchLengthsRecursively(TreeAln& traln, nodeptr p, double multiHere); 

  std::vector<Branch> storedBranches; 

} ; 
