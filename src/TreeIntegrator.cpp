#include "TreeIntegrator.hpp"
#include "proposals/BranchLengthMultiplier.hpp"
#include "GibbsProposal.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "Branch.hpp"
#include "priors/ExponentialPrior.hpp"
#include "proposals/TreeLengthMultiplier.hpp"
#include "ProposalRegistry.hpp"
#include "proposals/GibbsBranchLength.hpp"
#include "proposals/BranchLengthMultiplier.hpp"
#include "SprMove.hpp"
#include "Arithmetics.hpp"
#include "eval/FullCachePolicy.hpp"

// #define NUM_GEN_BASE 1000
// #define THINNING 10 

#define ESS_THRESH 100.



TreeIntegrator::TreeIntegrator(std::shared_ptr<TreeAln> tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed)
{
  auto && plcy = std::unique_ptr<ArrayPolicy>(new FullCachePolicy(*tralnPtr, true, true ));
  auto eval = LikelihoodEvaluator(*tralnPtr, plcy.get());

#ifdef DEBUG_LNL_VERIFY
  eval.setDebugTraln(debugTree);
#endif

  // s.t. we do not have to care about the branch length linking problem 
  assert(tralnPtr->getNumberOfPartitions() == 1 ); 


  std::vector<std::unique_ptr<AbstractParameter> > params; 
  params.emplace_back(std::unique_ptr<AbstractParameter>(new BranchLengthsParameter(0,0, {0}))); 
  for(nat i = 0; i < tralnPtr->getNumberOfPartitions(); ++i)
    params[0]->addPartition(i);

  double lambda = 10;
  params[0]->setPrior(std::unique_ptr< AbstractPrior>(new ExponentialPrior(lambda)));

  std::vector<std::unique_ptr<AbstractProposal> >  proposals;   

  // proposals.emplace_back(new BranchLengthMultiplier( ProposalRegistry::initBranchLengthMultiplier)); 
  proposals.emplace_back(new GibbsBranchLength()); 
  // proposals.emplace_back(new TreeLengthMultiplier(ProposalRegistry::initTreeLengthMultiplier)); 

  proposals[0]->addPrimaryParameter( std::move(std::unique_ptr<AbstractParameter>(params[0]->clone()))); 
  // proposals[1]->addPrimaryParameter( std::move(std::unique_ptr<AbstractParameter>(params[0]->clone()))); 
  // proposals[2]->addPrimaryParameter( std::move(std::unique_ptr<AbstractParameter>(params[0]->clone()))); 

  std::vector<ProposalSet> pSets; 

  integrationChain = std::unique_ptr<Chain>( new Chain(seed, tralnPtr, proposals, pSets, std::move(eval), false ));
  integrationChain->getEvaluator().imprint(*tralnPtr);  
}


void TreeIntegrator::prepareChain( const TreeAln &otherTree, bool copyEverything)
{
  auto &traln = integrationChain->getTraln();
  if(copyEverything)
    traln.copyModel(otherTree); 
  
  integrationChain->getEvaluator().evaluate(traln, traln.getAnyBranch(), true); 
  integrationChain->reinitPrior(); 
}


#include "TreeRandomizer.hpp"


// tree must be prepared! 
Branch2Stat TreeIntegrator::integrateTree(double essThresh, LikelihoodEvaluator &eval )
{
  // tout << "integrating tree " << numgen << " generations" << std::endl; 
  auto samples = std::unordered_map<BranchLength, std::vector<double>>{}; 
  auto params = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  auto& traln = *(integrationChain->getTralnPtr()); 
  assert(params.size( )== 1); 
  auto param = params[0]; 

  
  eval.evaluate(traln, traln.getAnyBranch(),true); 
  eval.freeMemory(); 
  integrationChain->setLikelihood(traln.getTrHandle().likelihood); 
  integrationChain->suspend(); 
  integrationChain->resume(true, false ); 

  bool converged = false; 
  while(not converged )
    {
      tout << "doing 10000 "  << std::endl; 
      
      for(nat i = 0; i < 10000; ++i)
	{
	  integrationChain->step();

	  if(i % ( traln.getNumberOfBranches( )* 5 )  == 0)
	    {
	      for(auto branch : traln.extractBranches(param))
		{
		  branch = traln.getBranch(branch.toPlain(), param); 
		  samples[branch].push_back(branch.getInterpretedLength(traln,param)); 
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
      // tout << "worst was: " << MAX_SCI_PRECISION << worstSample << std::endl; 
    }

  auto result = Branch2Stat{}; 
  for(auto sample : samples)  
    {
      result[sample.first] = make_pair(Arithmetics::getMean(sample.second), 
				       Arithmetics::getVariance(sample.second));
    }

  return result; 
}






void TreeIntegrator::integrateAllBranchesNew(const TreeAln &otherTree, std::string runid, nat sprDistance)
{
  auto&& ss = std::stringstream{}; 
  ss << PROGRAM_NAME ; 
  ss << "_sprIntegration-" << sprDistance << "." << runid; 
  auto&& outfile =  std::ofstream (ss.str());
  outfile << MAX_SCI_PRECISION ; 

  outfile
    << "moveId\t"
    << "move\t" 
    << "branch\t" 
    << "bl.mean.orig\t" 
    << "lnlRatio\t"  
    << "nniCat\t" 
    << "nniId\t"
    << "bl.ratio.nni\t"
    << "mrbCat\t"
    << "mrbId\t"
    << "bl.ratio.mrb" << std::endl; 
  
  prepareChain(otherTree, true); 

  auto& traln  = *(integrationChain->getTralnPtr()); 
  auto &eval = integrationChain->getEvaluator(); 
  eval.freeMemory();
  
  auto params = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  assert(params.size() == 1 );
  auto param = params[0]; 
 
  auto allBranchesBefore = traln.extractBranches();
  auto allBranchLengthsBefore = traln.extractBranches(param) ;

  auto initBranch2Stat = integrateTree(ESS_THRESH, eval); 
  double initLnl = traln.getTrHandle().likelihood;  

  auto moves = SprMove::getAllUniqueMoves(traln,sprDistance); 

  nat moveCtr = 0; 
  for(auto move : moves)
    {
      // tout << "integrating " << move << std::endl; 
      auto nameMapNni = move.getNames(traln, true); 
      auto nameMapMrb = move.getNames(traln, false); 

      move.applyToTree(traln, params);
      // eval.evaluate(traln, traln.getAnyBranch(), true); 
      // eval.freeMemory();


      auto branchStatAfter = integrateTree(ESS_THRESH, eval); 
      
      auto lnlAfter = traln.getTrHandle().likelihood; 

      nat blCtr = 0; 
      for(auto b : allBranchesBefore)
	{
	  auto beforeStat = initBranch2Stat.at(b.toBlDummy()); 
	  auto afterStatMrb = branchStatAfter.at(move.mapBranchSingleMapAfter(b).toBlDummy());
	  auto afterStatNni = branchStatAfter.at(move.mapBranchNniStepsAfter(b).toBlDummy()); 

	  auto nameTupleMrb = nameMapMrb[b]; 
	  auto nameTupleNni = nameMapNni[b]; 

	  outfile << moveCtr <<  "\t"
		  << move << "\t"
		  << b << "\t"
		  << beforeStat.first << "\t"
		  << lnlAfter - initLnl << "\t"
		  << nameTupleNni.first << "\t"
		  << nameTupleNni.second << "\t"
		  << afterStatNni.first / beforeStat.first << "\t"
		  << nameTupleMrb.first << "\t"
		  << nameTupleMrb.second << "\t"
		  << afterStatMrb.first / beforeStat.first << "\t"
		  << std::endl; 	  
	  ++blCtr; 
	}
      
      ++moveCtr; 
      move.revertTree(traln,params); 
    }
}
