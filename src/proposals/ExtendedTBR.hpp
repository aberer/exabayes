
#include "Path.hpp"
#include "AbstractProposal.hpp"
#include "TbrMove.hpp"

class ExtendedTBR : public AbstractProposal
{
public: 
  ExtendedTBR( double _extensionProb, double _multiplier); 
  virtual ~ExtendedTBR()  { }

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand); 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior); 
  virtual void resetState(TreeAln& traln, PriorBelief &prior); 
  virtual void autotune() {} 

  virtual AbstractProposal* clone() const; 

  virtual void readFromCheckpointCore(std::ifstream &in) {   } 
  virtual void writeToCheckpointCore(std::ofstream &out) { }  

private: 			// METHODS
  void drawPaths(TreeAln &traln, Randomness &rand); 
  void buildPath(Path &path, Branch bisectedBranch, TreeAln &traln, Randomness &rand );

private: 			// ATTRIBUTES
  double extensionProbability;   
  double multiplier; 
  TbrMove move; 
}; 
