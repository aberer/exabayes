/** 
    @file ParsimonySPR.hpp
    
    @brief implements a parsimony-biased SPR move similar to MrBayes. 

    @notice: MrBayes also reweights site patters -- currently we do
    not do that.
 */ 

#ifndef __PARSIMONY_SPR
#define __PARSIMONY_SPR

// why not evaluate stuff only on a few partitions and use that for guidance? 


#include <unordered_map>
#include "AbstractProposal.hpp"
#include "data-struct/Path.hpp"
#include "SprMove.hpp"
#include "eval/ParsimonyEvaluator.hpp"
#include "comm/RemoteComm.hpp"

typedef std::unordered_map<BranchPlain, double> weightMap; 
typedef std::unordered_map<BranchPlain,std::array<parsimonyNumber,2> > branch22states2score;


class ParsimonySPR : public AbstractProposal
{
public: 
  ParsimonySPR( double parsWarp, double blMulti, int depth); 
  virtual ~ParsimonySPR(){}

  virtual BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const; 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln) ; 
  virtual void autotune() ;

  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0) );}
  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out) const { } //disabled

  virtual std::vector<nat> getInvalidatedNodes(const TreeAln& traln) const; 

  virtual void printParams(std::ostream &out)  const ; 

  virtual AbstractProposal* clone() const; 

protected:			// METHODS
  weightMap getWeights(const TreeAln& traln,  branch22states2score insertions, LikelihoodEvaluator& eval) const; 
  branch22states2score determineScoresOfInsertions(TreeAln& traln, BranchPlain primeBranch, Randomness &rand, const branch22states2score &alreadyComputed ); 
  void traverse(const TreeAln &traln, nodeptr p, int distance ); 
  void testInsertParsimony(TreeAln &traln, nodeptr insertPos, nodeptr prunedTree, branch22states2score &result, int curDepth,  const branch22states2score& alreadyComputed); 

protected: 			// ATTRIBUTES
  double _parsWarp; 
  double _blMulti;   
  SprMove _move; 
  ParsimonyEvaluator _pEval;   
  BranchPlain _subtreeBranch; 
  int _depth;
  bool _parallelReduceAtEnd;

  static std::array<double,2> factors; 

  bool _proposePrimeBranch; 
}; 

#endif
