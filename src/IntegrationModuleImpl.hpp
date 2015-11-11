// just spliced this off for estetics. 

// experimental code that should not directly mix with production level code. 

#define STEPS_FOR_LNL 100
#define INTEGRATION_GENERATIONS 10000

// #define STEPS_FOR_LNL 10
// #define INTEGRATION_GENERATIONS 100
#define NR_STEPS 30


#include <sstream>
#include "priors/ExponentialPrior.hpp"
#include "proposals/BranchIntegrator.hpp"
#include "system/ProposalRegistry.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "eval/ParsimonyEvaluator.hpp"

#include "AdHocIntegrator.hpp"


void SampleMaster::branchLengthsIntegration(Randomness &rand)  
{
  assert(_runs.size() == 1 );   
  auto &run = _runs[0];   
  auto &chains = run.getChains(); 
  assert(chains.size() == 1); 
  auto &chain = chains[0]; 

  auto& traln  = chain.getTralnHandle(); 
  
  double lambda = 10 ; 
  
  auto&& ahInt = AdHocIntegrator (traln, nullptr ,rand.generateSeed(), _plPtr); 
  auto paramView = _runs[0].getChains()[0].getProposalView()[0]->getBranchLengthsParameterView();

  std::vector<std::pair<double,double> > parsAndMLBlen; 
  for(auto &branch : traln.extractBranches(ahInt.getBlParamView()[0]))
    {      
      tout << branch << std::endl; 
      auto initBranch = branch; 

      ahInt.prepareForBranch( branch.toPlain(), traln); 

      // sampling 
      auto samples = ahInt.integrate( branch.toPlain(), traln );

      double minHere = *( std::min_element(samples.begin(), samples.end()) ) ; 
      double maxHere = *( std::max_element(samples.begin(), samples.end()) ) ; 

      auto &&ss =stringstream{};
      ss << "samples." << _cl.getRunid()<< "." << branch.getPrimNode() << "-" << branch.getSecNode()   <<  ".tab" ;
      ofstream thisOut (ss.str());
      std::copy(samples.begin(), samples.end(), std::ostream_iterator<double>(thisOut, "\n")); 

      // lnl curve 
      ahInt.createLnlCurve(branch.toPlain(), _cl.getRunid(), traln, minHere, maxHere, STEPS_FOR_LNL);

      // optimization 
      double nrOpt = ahInt.printOptimizationProcess(branch,  _cl.getRunid(), lambda, NR_STEPS);

      // print parsimony length  
      double pLength = ahInt.getParsimonyLength(traln, branch.toPlain()); 
      parsAndMLBlen.push_back({pLength, nrOpt});

      // reset 
      traln.setBranch(initBranch, paramView[0]);
    }

  
  std::stringstream ss; 
  ss << "parsLengthVsML." << _cl.getRunid() << ".tab"; 
  std::ofstream out(ss.str());
  for(auto elem : parsAndMLBlen)
    out << elem.first << "\t" << elem.second << std::endl; 

  tout << "finished!" << endl; 
}
