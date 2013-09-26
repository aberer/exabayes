
#include "Path.hpp"
#include "AbstractProposal.hpp"
#include "TbrMove.hpp"

class ExtendedTBR : public AbstractProposal
{
public: 
  ExtendedTBR( double _extensionProb, double _multiplier); 
  virtual ~ExtendedTBR()  { }

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval); 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln); 
  virtual void resetState(TreeAln& traln); 
  virtual void autotune() {} 

  virtual AbstractProposal* clone() const; 

  // virtual Branch prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return Branch(0,0);}
  // virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}
    virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0)); }

  virtual void readFromCheckpointCore(std::istream &in) {   } 
  virtual void writeToCheckpointCore(std::ostream &out) const { }  

  virtual std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const; 

private: 			// METHODS
  void drawPaths(TreeAln &traln, Randomness &rand); 
  void buildPath(Path &path, BranchPlain bisectedBranch, TreeAln &traln, Randomness &rand );

private: 			// ATTRIBUTES
  double extensionProbability;   
  double multiplier; 
  TbrMove move; 
}; 
