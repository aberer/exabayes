#ifndef _LIKELIHOOD_SPR
#define _LIKELIHOOD_SPR


#include <unordered_map>

#include "AbstractProposal.hpp" 
#include "SprMove.hpp"
#include "model/Branch.hpp"

typedef std::unordered_map<BranchLength,double> BranchToLnlMap; 
typedef std::unordered_map<BranchLength, std::unordered_map<BranchLength, double>> BranchToNRD2; 


class LikelihoodSPR : public AbstractProposal
{
public: 
  LikelihoodSPR(nat minStep, nat maxStep, double likeWarp); 

  virtual BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const ; 
  
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln)  ; 
  virtual AbstractProposal* clone() const ;  
  virtual void autotune()   {}
    
  virtual void writeToCheckpointCore(std::ostream &out)const ;  
  virtual void readFromCheckpointCore(std::istream &in) ; 

  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand) { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0)); }
  
  virtual std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const; 

  void proposeBranches(TreeAln& traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval ); 
    
private: 			// METHODS 
  double scoreReattachment(TreeAln& traln, const BranchLength &reattachmentBranch, nodeptr prunedSubtree, LikelihoodEvaluator& eval, bool doRevert, bool isForward ) ; 
  BranchToLnlMap scoreOnBothSides(TreeAln &traln, const BranchLength &branch , nodeptr subtree, LikelihoodEvaluator &eval, bool isForward)  ; 
  void scoreReattachmentInRadius(TreeAln &traln, BranchLength attachment, nodeptr prunedSubtree, nat depth, LikelihoodEvaluator& eval, BranchToLnlMap &map, bool isForward)  ; 
  void determineMove(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval); 
  std::pair<BranchLength,double> drawReattachment(const BranchToLnlMap &map, Randomness &rand) const ; 
  BranchToLnlMap convertToProbMap(const BranchToLnlMap &map) const ; 

private: 
  SprMove move; 
  nat minStep ; 
  nat maxStep; 
  double likeWarp; 
  BranchToNRD2 map;
  BranchLength savedSubtreeBranch ; 
}; 


#endif
