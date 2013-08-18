#ifndef ADHOC_INTEGRATOR 
#define ADHOC_INTEGRATOR 

#include <sstream>
#include "priors/ExponentialPrior.hpp"
#include "proposals/BranchIntegrator.hpp"
#include "ProposalRegistry.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "ParsimonyEvaluator.hpp"
#include "RestoringLnlEvaluator.hpp"
#include "common.h"

// HACK 
struct noDeleter
{
  void operator ()(void *) {} 
};


class AdHocIntegrator
{
public: 
  AdHocIntegrator(std::shared_ptr<TreeAln> tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed); 
  std::vector<AbstractParameter*> getBlParamView() const ;
  void prepareForBranch( const Branch &branch, const TreeAln &traln); 
  std::vector<double> integrate(const Branch &branch, const TreeAln &otherTree, nat intGens, nat thinning); 
  void createLnlCurve(Branch branch, std::string runid, TreeAln & traln , double minHere, double maxHere, nat numSteps); 
  bool decideUponAcceptance(const TreeAln &traln, double prevLnl); 
  double getParsimonyLength(TreeAln &traln, const Branch &b );
  double printOptimizationProcess(const Branch& branch, std::string runid, double lambda, nat nrSteps); 
  void copyTree(const TreeAln &traln); 

private: 			// METHODS 


private: 			// ATTRIBUTES
  std::unique_ptr<Chain> integrationChain; 
}; 

#endif

