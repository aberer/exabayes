#ifndef _LIKELIHOOD_SPR 
#define _LIKELIHOOD_SPR 

#include "AbstractProposal.hpp"
#include "SprMove.hpp"

class LikelihoodSPR : public AbstractProposal
{
public: 
  explicit LikelihoodSPR(nat maxStep, double likeWarp);
  ~LikelihoodSPR(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln)  ; 
  virtual BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const ; 
  virtual std::vector<nat> getInvalidatedNodes(const TreeAln &traln) const ;  

  virtual void prepareForSetEvaluation( TreeAln &traln, LikelihoodEvaluator& eval) const  {} 
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0) ); }
  virtual AbstractProposal* clone() const {return new LikelihoodSPR(*this); }

  // NOOPS 
  virtual void writeToCheckpointCore(std::ostream &out)const {}
  virtual void readFromCheckpointCore(std::istream &in) {}
  virtual void autotune() {}

  virtual void printParams(std::ostream &out)  const ; 

private: 			// METHODS
    auto computeLikelihoodsOfInsertions(TreeAln &traln, LikelihoodEvaluator &eval, const BranchPlain& prunedTree, std::vector<AbstractParameter*> &params)  
      -> std::unordered_map<BranchPlain, log_double>; 
  auto  determineSprMove(TreeAln &traln, LikelihoodEvaluator &eval, Randomness &rand, const BranchPlain& prunedTree, std::vector<AbstractParameter*> &params)  
    -> std::tuple<SprMove, double>; 
  auto transformLikelihoods(std::unordered_map< BranchPlain, log_double > branch2lnl ) const 
    -> std::unordered_map<BranchPlain,double> ; 
  

private: 			// ATTRIBUTES
  nat _maxStep; 
  double _likeWarp; 
  SprMove _move; 
}; 

#endif

