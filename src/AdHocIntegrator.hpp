#ifndef ADHOC_INTEGRATOR 
#define ADHOC_INTEGRATOR 

#include <sstream>
#include "priors/ExponentialPrior.hpp"
#include "proposals/BranchIntegrator.hpp"
#include "ProposalRegistry.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "ParsimonyEvaluator.hpp"
#include "RestoringLnlEvaluator.hpp"


struct noDeleter
{
  void operator ()(void *) {} 
};


class AdHocIntegrator
{
public: 
  AdHocIntegrator(std::shared_ptr<TreeAln>  tralnPtr, randCtr_t seed ); 
  std::vector<AbstractParameter*> getBlParamView() const ;
  void prepareForBranch( Branch branch, TreeAln &traln); 
  std::pair<double,double> integrate(TreeAln &traln, std::string runid, Branch branch, nat intGens); 
  void createLnlCurve(Branch branch, std::string runid, TreeAln & traln , double minHere, double maxHere, nat numSteps); 
  double getParsimonyLength(TreeAln &traln, const Branch &b );
  double printOptimizationProcess(Branch branch, TreeAln &traln, std::string runid, double lambda, nat nrSteps); 

  void copyTree(const TreeAln &traln); 

private: 			// METHODS 
  std::pair<double,double> getMeanAndVar (const std::vector<double> &data ); 

private: 			// ATTRIBUTES
  std::unique_ptr<Chain> integrationChain; 
}; 

#endif

