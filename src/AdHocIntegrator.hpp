#ifndef ADHOC_INTEGRATOR 
#define ADHOC_INTEGRATOR 

#include <sstream>
#include "priors/ExponentialPrior.hpp"
#include "proposals/BranchIntegrator.hpp"
#include "system/ProposalRegistry.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "eval/ParsimonyEvaluator.hpp"
#include "common.h"

class Communicator; 

// HACK 
struct noDeleter
{
  void operator ()(void *) {} 
};


class AdHocIntegrator
{
public: 
  AdHocIntegrator(TreeAln& tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed, ParallelSetup* pl); 
  double getParsimonyLength(TreeAln &traln, const BranchPlain &b, Communicator& comm ); 
  std::vector<AbstractParameter*> getBlParamView() const ;
  void prepareForBranch( const BranchPlain &branch, const TreeAln &traln); 
  std::vector<double> integrate(const BranchPlain &branch, const TreeAln &otherTree); 
  void createLnlCurve(BranchPlain branch, std::string runid, TreeAln & traln , double minHere, double maxHere, nat numSteps); 
  double printOptimizationProcess(const BranchLength& branch, std::string runid, nat nrSteps, Communicator& comm); 
  void copyTree(const TreeAln &traln); 

private: 			// METHODS 


private: 			// ATTRIBUTES
  std::unique_ptr<Chain> integrationChain; 
}; 

#endif

