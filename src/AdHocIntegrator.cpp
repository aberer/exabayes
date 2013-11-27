#include "AdHocIntegrator.hpp"
#include "eval/ArrayReservoir.hpp"
#include "eval/ArrayRestorer.hpp"
#include "Arithmetics.hpp"
#include "eval/FullCachePolicy.hpp"


AdHocIntegrator::AdHocIntegrator(TreeAln &traln, std::shared_ptr<TreeAln> debugTree, randCtr_t seed)
{
  auto && plcy = std::unique_ptr<ArrayPolicy>(new FullCachePolicy(traln, true, true));
  auto&& res = std::make_shared<ArrayReservoir>(false);
  auto eval = LikelihoodEvaluator(traln, plcy.get() , res);

#ifdef DEBUG_LNL_VERIFY
  eval.setDebugTraln(debugTree);
#endif
  
  // s.t. we do not have to care about the branch length linking problem 
  assert(traln.getNumberOfPartitions() == 1 ); 

  std::vector<std::unique_ptr<AbstractParameter> > params; 
  params.emplace_back(std::unique_ptr<AbstractParameter>(new BranchLengthsParameter(0,0, {0}))); 
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    params[0]->addPartition(i);

  double lambda = 10;

  params[0]->setPrior(std::unique_ptr<AbstractPrior>(new ExponentialPrior(lambda)));

  vector<unique_ptr<AbstractProposal> >  proposals;   
  proposals.emplace_back( new BranchIntegrator (ProposalRegistry::initBranchLengthMultiplier)); 
  proposals[0]->addPrimaryParameter( std::move(params[0])); 

  std::vector<ProposalSet> pSets; 

  integrationChain = std::unique_ptr<Chain>( new Chain(seed, traln, proposals, pSets, std::move(eval), false ));
  integrationChain->getEvaluator().imprint(traln);
}


void AdHocIntegrator::copyTree(const TreeAln &traln)
{
  auto &myTree = integrationChain->getTralnHandle();

  myTree = traln; 
  
  auto& eval = integrationChain->getEvaluator();
  eval.evaluate(myTree, traln.getAnyBranch(), true , true);
}


void AdHocIntegrator::prepareForBranch( const BranchPlain &branch,  const TreeAln &otherTree)
{
  copyTree(otherTree); 

  auto &traln = integrationChain->getTralnHandle();

  auto ps = integrationChain->getProposalView(); 
  auto paramView = ps[0]->getBranchLengthsParameterView();
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0]); 
  integrator->setToPropose(branch);      

  // std::cout << "important TODO: has the branch been reset? will not do that here" << std::endl; 

  integrationChain->getEvaluator().evaluate(traln, branch, true, true); 
  integrationChain->reinitPrior(); 
}


std::vector<AbstractParameter*> AdHocIntegrator::getBlParamView() const
{
  auto ps = integrationChain->getProposalView(); 
  auto paramView = ps[0]->getBranchLengthsParameterView();
  assert(ps.size() == 1 );   
  return paramView; 
}


std::vector<double> AdHocIntegrator::integrate( const BranchPlain &branch, const TreeAln &otherTree)
{
  auto result =  std::vector<double>{}; 
  auto& traln = integrationChain->getTralnHandle(); 

  // tout << "integrating " << branch << std::endl; 

  prepareForBranch(branch, otherTree); 
  auto paramView = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  assert(paramView.size( )== 1 ); 

  auto backup = otherTree.getBranch(branch, paramView[0]); 

  bool converged = false; 
  while(not converged)
    {
      for(nat i = 0; i < 10000; ++i) 
	{	  
	  integrationChain->step();
	  auto elem = traln.getBranch(branch, paramView[0]); 
	  auto iLen = elem.getInterpretedLength(traln, paramView[0]);
	  if (i % 10 == 0)
	    result.push_back(iLen); ; 
	}
      
      auto ess = Arithmetics::getEffectiveSamplingSize(result); 
      converged = ess > 200; 
    }

  traln.setBranch(backup, paramView[0]); 
  return result;   
}


/** 
    @brief gets the optimimum 
 */ 
double AdHocIntegrator::printOptimizationProcess(const BranchLength& branch, std::string runid, double lambda, nat nrSteps)
{
  auto &traln = integrationChain->getTralnHandle(); 

  auto paramView = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  auto tmpBranch = branch; 
  // tout << "optimizing the branch using nr" << endl; 
  std::stringstream ss; 
  ss << "nr-length." << runid  << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
  std::ofstream thisOut(ss.str()); 

  double result = 0;  
  double curVal = 0.1; 
  tmpBranch.setConvertedInternalLength(traln, paramView[0], curVal); 
  double secDerivative = 0; 
  double firstDerivative = 0; 

  double prevVal = curVal; 
  curVal = tmpBranch.getLength();
  
  for(nat i = 0; i < nrSteps; ++i )
    {
#if HAVE_PLL != 0 
      makenewzGeneric(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &curVal, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE, NULL); 
#else 
      makenewzGeneric(&(traln.getTrHandle()), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &curVal, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE, NULL); 	  
#endif
      tmpBranch.setLength(result);
      thisOut << prevVal <<  "\t" << firstDerivative << "\t" << secDerivative << endl; 	
      prevVal = tmpBranch.getInterpretedLength(traln, paramView[0]); 
      curVal = result; 
    } 

  tmpBranch.setConvertedInternalLength(traln, paramView[0], prevVal); 
  double something = tmpBranch.getLength(); 

#if HAVE_PLL != 0
  makenewzGeneric(&traln.getTrHandle(), &traln.getPartitionsHandle(), 
		  branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		  &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE, NULL); 
#else 
  makenewzGeneric(&traln.getTrHandle(), 
		  branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		  &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE, NULL); 
#endif
      
  // thisOut << prevVal << "\t" << firstDerivative << "\t" << secDerivative << endl; 

  thisOut.close(); 

  return prevVal; 
}



void AdHocIntegrator::createLnlCurve(BranchPlain branch, std::string runid, TreeAln & traln , double minHere, double maxHere, nat numSteps)
{
  auto paramView  = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  assert(paramView.size( )== 1 ); 
  auto param = paramView[0]; 
  std::stringstream ss; 
  ss << "lnl." << runid << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
  std::ofstream thisOut(ss.str());

  auto& eval = integrationChain->getEvaluator();

  if(maxHere != minHere)
    {
      for(double i = minHere; i < maxHere+0.00000001 ; i+= (maxHere-minHere)/ numSteps)
	{
	  auto b = traln.getBranch(branch, param); 
	  b.setConvertedInternalLength(traln, paramView[0], i); 
	  traln.setBranch(b, paramView[0]); 
	  eval.evaluate(traln, branch, false, true);
	  double lnl = traln.getTrHandle().likelihood; 
	  thisOut << i << "\t" << setprecision(std::numeric_limits<double>::digits10) << lnl << endl; 
	}
    }
  else
    thisOut << minHere << "\t" << "NA" << endl; 
  thisOut.close(); 
}


double AdHocIntegrator::getParsimonyLength(TreeAln &traln, const BranchPlain &b )
{
  auto pEval =  ParsimonyEvaluator{}; 
  auto branchLength = std::vector<nat>{} ; 
  // auto state2pars = 
  pEval.evaluate(traln, b.findNodePtr(traln), true, true  );

  assert(traln.getNumberOfPartitions() == 1 ); 
  auto& partition =  traln.getPartition(0);
  auto length = partition.upper - partition.lower; 
  
  // TODO this is incorrect!
  assert(0); 

  double result =  double(branchLength[0])  / double(length); 
  return result; 
}

