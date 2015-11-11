#include "TreeIntegrator.hpp"
#include "extensions.hpp" 
#include "ArrayReservoir.hpp"
#include "BranchLengthsParameter.hpp"
#include "ChosenBranchIntegrator.hpp"
#include "Branch.hpp"
#include "ExponentialPrior.hpp"
#include "TreeLengthMultiplier.hpp"
#include "ProposalRegistry.hpp"

#include "MoveEnumeration.hpp"

#include "SprMove.hpp"
#include "Arithmetics.hpp"
#include "FullCachePolicy.hpp"

#include "TreeRandomizer.hpp"

#include "BranchBackup.hpp"

// #define NUM_GEN_BASE 1000
// #define THINNING 10 

#define ESS_THRESH 200
#define MAX_MOVES 200  

using std::get; 
using std::endl; 

TreeIntegrator::TreeIntegrator(TreeAln& traln, std::shared_ptr<TreeAln> debugTree, randCtr_t seed, ParallelSetup* plPtr)
  : integrationChain{nullptr}
{
  auto && plcy = make_unique<FullCachePolicy>(traln, true, true );
  auto res = std::make_shared<ArrayReservoir>(false);
  auto eval = LikelihoodEvaluator(traln, plcy.get(), res, plPtr);

#ifdef DEBUG_LNL_VERIFY
  eval.setDebugTraln(debugTree);
#endif

  // s.t. we do not have to care about the branch length linking problem 
  assert(traln.getNumberOfPartitions() == 1 ); 


  auto params = std::vector<std::unique_ptr<AbstractParameter> >{} ; 
  params.emplace_back(std::unique_ptr<AbstractParameter>(new BranchLengthsParameter(0,0, {0}))); 
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    params[0]->addPartition(i);

  double lambda = 10;
  params[0]->setPrior(std::unique_ptr< AbstractPrior>(new ExponentialPrior(lambda)));

  auto &&proposals  = std::vector<std::unique_ptr<AbstractProposal> > {};   

  proposals.emplace_back( new ChosenBranchIntegrator(ProposalRegistry::initBranchLengthMultiplier));

  proposals[0]->addPrimaryParameter( params[0]->getId() ); 

  auto &&paramList = ParameterList(std::move(params)); 
  integrationChain = make_unique<Chain>( seed, traln, std::move(paramList), proposals, std::vector<ProposalSet>{}, std::move(eval), false );
  integrationChain->getEvaluator().imprint(traln);  
}


void TreeIntegrator::prepareChain( const TreeAln &otherTree, bool copyEverything)
{
  auto &traln = integrationChain->getTralnHandle();
  if(copyEverything)
    traln = otherTree; 

  integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(), true); 
  integrationChain->reinitPrior(); 
}





// tree must be prepared! 
std::vector<BranchStat> TreeIntegrator::integrateTree(double essThresh, LikelihoodEvaluator &eval, std::vector<BranchPlain> branches)
{
  auto samples = std::unordered_map<BranchPlain, std::vector<double>>{}; 
  auto params = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  auto& traln = integrationChain->getTralnHandle(); 
  assert(params.size( )== 1); 
  auto param = params[0]; 

  eval.evaluate(traln, traln.getAnyBranch(),true); 
  eval.freeMemory(); 
  integrationChain->setLikelihood(  traln.getLikelihood()); 
  integrationChain->suspend(); 
  integrationChain->resume( ); 
  
  auto proposals = integrationChain->getProposalView();
  assert(proposals.size() == 1); 
  auto proposal = proposals[0]; 
  assert(dynamic_cast<ChosenBranchIntegrator*>(proposal) != nullptr);
  auto casted =  dynamic_cast<ChosenBranchIntegrator*>(proposal); 
  casted->setChosenBranches(branches); 

  bool converged = false; 
  while(not converged )
    {
      for(nat i = 0; i < 10000; ++i)
	{
	  integrationChain->step();

	  if(i % ( branches.size() * 10 )  == 0)
	    {
	      for(auto branch : branches)
		{
		  auto b  = traln.getBranch(branch, param); 
		  samples[b.toPlain()].push_back(b.getInterpretedLength(param)); 
		}
	    }
	}

      converged = true; 
      
      double least = std::numeric_limits<double>::infinity(); 
      double mean = 0.; 

      auto worstSample = std::vector<double> {}; 
      for(auto &elem : samples)
	{
	  double ess = Arithmetics::getEffectiveSamplingSize(elem.second); 

	  if(std::isnan(ess))
	    {
	      tout << "for elem "  << elem.first << " the ess was " << ess << "\tsamples: " << elem.second << std::endl; 
	    }
	  
	  converged &= ess > ESS_THRESH ; 
	  mean += ess; 
	  if(ess < least)
	    {
	      least = ess; 
	      worstSample = elem.second;  
	    }
	}
      
      tout << "least ESS: " << SOME_FIXED_PRECISION << least << "\tmean(ess)=" << mean / double(traln.getNumberOfBranches())  << std::endl; 
    }

  auto result = std::vector<BranchStat>{};
  for(auto sample : samples)  
    result.emplace_back(sample.first, Arithmetics::getMean(sample.second), Arithmetics::getVariance(sample.second)); 

  return result; 
}


void TreeIntegrator::integrateAllMoves(const TreeAln &otherTree, std::string runid, nat sprDistance)
{
  std::srand ( unsigned ( std::time(0) ) );

  // get file names 
  auto &&posteriorFile = std::ofstream{std::string{PROGRAM_NAME } + "_posteriors." + runid}; 
  posteriorFile << MAX_SCI_PRECISION ; 
  auto &&moveFile = std::ofstream{std::string{PROGRAM_NAME} + "_moves." + runid}; 
  moveFile << MAX_SCI_PRECISION ; 
  // end 

  prepareChain(otherTree, true); 

  auto& traln  = integrationChain->getTralnHandle(); 
  auto &eval = integrationChain->getEvaluator(); 
  eval.freeMemory();
  
  auto params = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  assert(params.size() == 1 );
  auto param = params[0]; 

  // save init state 
  auto allBranchesBefore = traln.extractBranches(param);
  auto initBackup = BranchBackup{};
  for(auto &v : allBranchesBefore)
    initBackup.extend(traln, param, v);
  
  eval.evaluate(traln,traln.extractBranches()[0], true);
  auto initLnl = traln.getLikelihood();

  moveFile << "INIT\t"  << "NA\t" << initLnl << std::endl; 

  tout << "initial integration"<< std::endl; 
  auto initBranch2Stat = integrateTree(ESS_THRESH, eval,traln.extractBranches() ); 
  for(auto elem : initBranch2Stat )
    {
      auto b =  BranchPlain(0,0); 
      auto mean = 0.; 
      auto var = 0.; 
      std::tie(b,mean,var) = elem; 
      posteriorFile << "INIT\t" << b.getPrimNode() << "\t" << b.getSecNode() << "\t" << mean << "\t" << std::sqrt(var) << std::endl; 
    }
  tout << "DONE!" << std::endl; 

  nat moveCtr = 0; 

  auto dist = 0; 
  auto res = std::string(getEnvironmentVariable("MOVE_DIST")); 
  dist = atoi(res.c_str());
  assert(dist != 0); 
  
  auto allMoves = MoveEnumeration::getAllUniqueMoves(traln,dist); 

  if(allMoves.size() > MAX_MOVES )
    {
      tout << "using only " << MAX_MOVES <<  " out of " << allMoves.size() << " moves of length "  << sprDistance << std::endl; 
      std::random_shuffle(begin(allMoves), end(allMoves) );
      allMoves = std::vector<SprMove>(begin(allMoves) , begin(allMoves) + MAX_MOVES ); 
    }
      
      
  for(auto move : allMoves)
    {
      initBackup.resetFromBackup(traln);
      integrationChain->reinitPrior();
      move.apply(traln, params);

      tout << "integrating " << move << std::endl; 

      auto relevantBranches = move.getInverseMove().getBranchesByDistance(traln,4); 
      auto relBraList =  std::vector<BranchPlain>(begin(relevantBranches), end(relevantBranches)); 

      tout << "relevant branches are " << relBraList << std::endl; 
      
      auto branchStatAfter = integrateTree(ESS_THRESH, eval, relBraList); 

      eval.evaluate(traln,relBraList.at(0), true); 
      auto afterLnl = traln.getLikelihood();	
      moveFile << moveCtr << "\t" << move << "\t"  << afterLnl << std::endl ;

      for(auto elem : branchStatAfter)
	{
	  auto b = BranchPlain(0,0); 
	  auto mean = 0.; 
	  auto var = 0.; 
	  std::tie(b,mean,var) = elem ; 
	  posteriorFile << moveCtr << "\t"  << b.getPrimNode() << "\t" << b.getSecNode() << "\t" << mean << "\t" << std::sqrt(var) << std::endl; 
	}

      ++moveCtr; 
      move.getInverseMove().apply(traln,params); 
    }
}
