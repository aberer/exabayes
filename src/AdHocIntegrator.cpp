#include "AdHocIntegrator.hpp"
#include "ArrayRestorer.hpp"



AdHocIntegrator::AdHocIntegrator(std::shared_ptr<TreeAln>  tralnPtr, randCtr_t seed)
{
  auto restorer = make_shared<ArrayRestorer>(*tralnPtr);  
  auto eval = std::unique_ptr<LikelihoodEvaluator>( new RestoringLnlEvaluator(restorer)); 
  
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

  integrationChain = std::unique_ptr<Chain>( new Chain(seed, tralnPtr, proposals, pSets, std::move(eval) ));
  integrationChain->getEvaluator()->imprint(*tralnPtr);
}



std::pair<double,double> AdHocIntegrator::getMeanAndVar (const std::vector<double> &data )
{
  double mean = 0; 
  for(auto d: data)
    mean += d; 
  mean /= data.size(); 
  
  double var = 0; 
  for(auto d : data)
    var += pow(d - mean, 2); 
  var /= data.size(); 

  return std::pair<double,double>(mean,var);  
}


void AdHocIntegrator::copyTree(const TreeAln &traln)
{
  TreeAln &myTree = integrationChain->getTraln();
  myTree.copyModel(traln);
  auto eval = integrationChain->getEvaluator();
  auto tr = myTree.getTr(); 

  eval->evaluate(myTree, Branch(tr->start->number,tr->start->back->number ), true );

  // tout << "copied tree " << traln << std::endl; 
  tout << "copied tree " << TreePrinter(true, true ,false).printTree(myTree) << std::endl; 
  // exit(0); 
}


void AdHocIntegrator::prepareForBranch( Branch branch, TreeAln &traln)
{
  copyTree(traln); 

  auto ps = integrationChain->getProposalView(); 
  auto paramView = ps[0]->getBranchLengthsParameterView();
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0]); 
  integrator->setToPropose(branch);      

  tout << "my param view is " << paramView << std::endl;

  // setup 
  traln.setBranch(branch, paramView); 

  integrationChain->getEvaluator()->evaluate(traln, branch, true); 
  integrationChain->reinitPrior(); 
}


std::vector<AbstractParameter*> AdHocIntegrator::getBlParamView() const
{
  auto ps = integrationChain->getProposalView(); 
  auto paramView = ps[0]->getBranchLengthsParameterView();
  assert(ps.size() == 1 );   
  return paramView; 
}


std::pair<double,double> AdHocIntegrator::integrate(TreeAln &traln, std::string runid, Branch branch, nat intGens)
{
  double minHere = 1000,
    maxHere = 0; 

  auto paramView = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();

  ////////////////////////
  // integrate branch   //
  ////////////////////////      
  stringstream ss; 
  ss << "samples." << runid<< "." << branch.getPrimNode() << "-" << branch.getSecNode()   <<  ".tab" ;
  ofstream thisOut (ss.str());       

  // run the chain to integrate 
  // tout << "integrating branch " << branch << endl; 
  for(nat i = 0; i < intGens; ++i) 
    {	  
      integrationChain->step();
      auto elem = traln.getBranch(branch.findNodePtr(traln), paramView); 
      auto iLen = elem.getInterpretedLength(traln, paramView[0]);
      if (i % 10 == 0)
	thisOut << iLen << endl; 
      if(iLen < minHere)
	minHere = iLen; 
      if(maxHere < iLen)
	maxHere = iLen; 
    }      

  thisOut.close();
  
  return std::pair<double,double>(minHere, maxHere);    
}


/** 
    @brief gets the optimimum 
 */ 
double AdHocIntegrator::printOptimizationProcess(Branch branch, TreeAln &traln, std::string runid, double lambda, nat nrSteps)
{
  auto paramView = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  Branch tmpBranch = branch; 
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
  curVal = tmpBranch.getLength(paramView[0]);
  
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
      tmpBranch.setLength(result, paramView[0]);
      thisOut << prevVal <<  "\t" << firstDerivative << "\t" << secDerivative << endl; 	
      prevVal = tmpBranch.getInterpretedLength(traln, paramView[0]); 
      curVal = result; 
    } 

  tmpBranch.setConvertedInternalLength(traln, paramView[0], prevVal); 
  double something = tmpBranch.getLength(paramView[0]); 

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



void AdHocIntegrator::createLnlCurve(Branch branch, std::string runid, TreeAln & traln , double minHere, double maxHere, nat numSteps)
{
  auto paramView  = integrationChain->getProposalView()[0]->getBranchLengthsParameterView();
  std::stringstream ss; 
  ss << "lnl." << runid << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
  std::ofstream thisOut(ss.str());

  auto eval = integrationChain->getEvaluator();

  if(maxHere != minHere)
    {
      for(double i = minHere; i < maxHere+0.00000001 ; i+= (maxHere-minHere)/ numSteps)
	{
	  Branch b = branch; 
	  b.setConvertedInternalLength(traln, paramView[0], i); 
	  traln.setBranch(b, paramView[0]); 
	  eval->evaluate(traln, branch, false);
	  double lnl = traln.getTr()->likelihood; 
	  thisOut << i << "\t" << setprecision(std::numeric_limits<double>::digits10) << lnl << endl; 
	}
    }
  else
    thisOut << minHere << "\t" << "NA" << endl; 
  thisOut.close(); 
}


double AdHocIntegrator::getParsimonyLength(TreeAln &traln, const Branch &b )
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

