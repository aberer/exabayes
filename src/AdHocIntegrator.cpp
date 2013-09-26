#include "AdHocIntegrator.hpp"
#include "eval/ArrayRestorer.hpp"
#include "Arithmetics.hpp"
#include "eval/FullCachePolicy.hpp"


AdHocIntegrator::AdHocIntegrator(std::shared_ptr<TreeAln>  tralnPtr, std::shared_ptr<TreeAln> debugTree, randCtr_t seed)
{
  // auto eval = std::unique_ptr<LikelihoodEvaluator>( new RestoringLnlEvaluator(*tralnPtr)); 
  auto && plcy = std::unique_ptr<ArrayPolicy>(new FullCachePolicy(*tralnPtr, true, true));
  auto eval = LikelihoodEvaluator(*tralnPtr, plcy.get() );

#ifdef DEBUG_LNL_VERIFY
  eval.setDebugTraln(debugTree);
#endif
  
  // s.t. we do not have to care about the branch length linking problem 
  assert(tralnPtr->getNumberOfPartitions() == 1 ); 

  std::vector<std::unique_ptr<AbstractParameter> > params; 
  params.emplace_back(std::unique_ptr<AbstractParameter>(new BranchLengthsParameter(0,0))); 
  for(nat i = 0; i < tralnPtr->getNumberOfPartitions(); ++i)
    params[0]->addPartition(i);

  double lambda = 10;

  params[0]->setPrior(make_shared< ExponentialPrior>(lambda));

  vector<unique_ptr<AbstractProposal> >  proposals;   
  proposals.emplace_back( new BranchIntegrator (ProposalRegistry::initBranchLengthMultiplier)); 
  proposals[0]->addPrimaryParameter( std::move(params[0])); 

  std::vector<ProposalSet> pSets; 

  integrationChain = std::unique_ptr<Chain>( new Chain(seed, tralnPtr, proposals, pSets, std::move(eval), false ));
  integrationChain->getEvaluator().imprint(*tralnPtr);
}


void AdHocIntegrator::copyTree(const TreeAln &traln)
{
  TreeAln &myTree = integrationChain->getTraln();
  myTree.copyModel(traln);
  auto& eval = integrationChain->getEvaluator();
  auto tr = myTree.getTr(); 
  eval.evaluate(myTree, BranchPlain(tr->start->number,tr->start->back->number ), true );
}


void AdHocIntegrator::prepareForBranch( const BranchPlain &branch,  const TreeAln &otherTree)
{
  copyTree(otherTree); 

  auto &traln = integrationChain->getTraln();

  auto ps = integrationChain->getProposalView(); 
  auto paramView = ps[0]->getBranchLengthsParameterView();
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0]); 
  integrator->setToPropose(branch);      

  std::cout << "important TODO: has the branch been reset? will not do that here" << std::endl; 

  // setup 
  // traln.setBranch(branch, paramView); 

  integrationChain->getEvaluator().evaluate(traln, branch, true); 
  integrationChain->reinitPrior(); 
}


std::vector<AbstractParameter*> AdHocIntegrator::getBlParamView() const
{
  auto ps = integrationChain->getProposalView(); 
  auto paramView = ps[0]->getBranchLengthsParameterView();
  assert(ps.size() == 1 );   
  return paramView; 
}


bool AdHocIntegrator::decideUponAcceptance(const TreeAln &traln, double prevLnl)
{
  copyTree(traln); 
  double nowLnl = integrationChain->getTraln().getTr()->likelihood; 
  double acc = integrationChain->getChainRand().drawRandDouble01(); 
  
  return acc < exp(nowLnl - prevLnl); 
}




std::vector<double> AdHocIntegrator::integrate( const BranchPlain &branch, const TreeAln &otherTree , nat intGens, nat thinning)
{
  auto result =  std::vector<double>{}; 
  auto& traln = integrationChain->getTraln(); 

  // tout << "integrating " << branch << std::endl; 

  prepareForBranch(branch, otherTree); 
  auto paramView = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  assert(paramView.size( )== 1 ); 

  auto backup = otherTree.getBranch(branch, paramView[0]); 
  for(nat i = 0; i < intGens; ++i) 
    {	  
      integrationChain->step();
      auto elem = traln.getBranch(branch, paramView[0]); 
      auto iLen = elem.getInterpretedLength(traln, paramView[0]);
      if (i % thinning == 0)
	result.push_back(iLen); ; 
    }      

  traln.setBranch(backup, paramView[0]); 

  return result;   
}


/** 
    @brief gets the optimimum 
 */ 
double AdHocIntegrator::printOptimizationProcess(const BranchLength& branch, std::string runid, double lambda, nat nrSteps)
{
  auto &traln = integrationChain->getTraln(); 

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
      makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &curVal, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
#else 
      makenewzGeneric(traln.getTr(), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &curVal, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 	  
#endif
      tmpBranch.setLength(result);
      thisOut << prevVal <<  "\t" << firstDerivative << "\t" << secDerivative << endl; 	
      prevVal = tmpBranch.getInterpretedLength(traln, paramView[0]); 
      curVal = result; 
    } 

  tmpBranch.setConvertedInternalLength(traln, paramView[0], prevVal); 
  double something = tmpBranch.getLength(); 

#if HAVE_PLL != 0
  makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), 
		  branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		  &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
#else 
  makenewzGeneric(traln.getTr(), 
		  branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		  &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
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
	  eval.evaluate(traln, branch, false);
	  double lnl = traln.getTr()->likelihood; 
	  thisOut << i << "\t" << setprecision(std::numeric_limits<double>::digits10) << lnl << endl; 
	}
    }
  else
    thisOut << minHere << "\t" << "NA" << endl; 
  thisOut.close(); 
}


double AdHocIntegrator::getParsimonyLength(TreeAln &traln, const BranchPlain &b )
{
  ParsimonyEvaluator pEval; 
  std::vector<nat> partitionParsimony; 
  std::vector<nat> branchLength; 
  pEval.evaluate(traln, b.findNodePtr(traln), true, partitionParsimony, branchLength );

  assert(traln.getNumberOfPartitions() == 1 ); 
  auto *partition =  traln.getPartition(0);
  auto length = partition->upper - partition->lower; 

  double result =  double(branchLength[0])  / double(length); 
  return result; 
}

