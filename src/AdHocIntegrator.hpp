#ifndef ADHOC_INTEGRATOR 
#define ADHOC_INTEGRATOR 

#include <sstream>
#include "priors/ExponentialPrior.hpp"
#include "proposals/BranchIntegrator.hpp"
#include "ProposalRegistry.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "ParsimonyEvaluator.hpp"
#include "common.h"

// HACK 
struct noDeleter
{
  void operator ()(void *) {} 
};


class AdHocIntegrator
{
public: 
  AdHocIntegrator(TreeAln& tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed); 
  std::vector<AbstractParameter*> getBlParamView() const ;
  void prepareForBranch( const BranchPlain &branch, const TreeAln &traln); 
  std::vector<double> integrate(const BranchPlain &branch, const TreeAln &otherTree); 
  void createLnlCurve(BranchPlain branch, std::string runid, TreeAln & traln , double minHere, double maxHere, nat numSteps); 
  // bool decideUponAcceptance(const TreeAln &traln, double prevLnl); 
  double getParsimonyLength(TreeAln &traln, const BranchPlain &b );
  double printOptimizationProcess(const BranchLength& branch, std::string runid, double lambda, nat nrSteps); 
  void copyTree(const TreeAln &traln); 

private: 			// METHODS 


private: 			// ATTRIBUTES
  std::unique_ptr<Chain> integrationChain; 
}; 

#endif

