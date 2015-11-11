
#include "data-struct/Path.hpp"
#include "AbstractProposal.hpp"
#include "TbrMove.hpp"

class ExtendedTBR : public AbstractProposal
{
public: 
  ExtendedTBR( double _extensionProb, double _multiplier); 

  BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const ; 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval); 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion); 
  virtual void resetState(TreeAln& traln); 
  virtual void autotune() {} 

  virtual AbstractProposal* clone() const; 

  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0)); }

  virtual void readFromCheckpointCore(std::istream &in) {   } 
  virtual void writeToCheckpointCore(std::ostream &out) const { }  

  virtual std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const; 

private: 			// METHODS
  void drawPaths(TreeAln &traln, Randomness &rand); 
  void buildPath(Path &path, BranchPlain bisectedBranch, TreeAln &traln, Randomness &rand );

private: 			// ATTRIBUTES
  double _extensionProbability;   
  double _multiplier; 
  TbrMove _move; 
}; 
