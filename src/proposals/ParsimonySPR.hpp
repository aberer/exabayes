/** 
    @file ParsimonySPR.hpp
    
    @brief implements a parsimony-biased SPR move similar to MrBayes. 

    @notice: MrBayes also reweights site patters -- currently we do
    not do that.
 */ 

#ifndef __PARSIMONY_SPR
#define __PARSIMONY_SPR

#include <unordered_map>

#include "axml.h"
#include "AbstractProposal.hpp"
#include "Path.hpp"
#include "SprMove.hpp"

typedef unordered_map<Branch, double, BranchHashNoLength, BranchEqualNoLength> weightMap; 
typedef unordered_map<Branch,vector<nat>, BranchHashNoLength, BranchEqualNoLength> scoreMap; 

class ParsimonySPR : public AbstractProposal
{
public: 
  ParsimonySPR( double parsWarp, double blMulti); 
  virtual ~ParsimonySPR(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune() ;

  AbstractProposal* clone() const; 

protected: 
  double parsWarp; 
  double blMulti;   
  Path path; 

  weightMap getWeights(const TreeAln& traln, const scoreMap &insertions) const; 
  // weightMap getForwardWeights( const scoreMap& insertions) const; 
  void determineSprPath(TreeAln& traln, Randomness &rand, double &hastings, PriorBelief &prior ); 
  
  SprMove move; 
  
}; 


#endif