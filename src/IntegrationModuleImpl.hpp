// just spliced this off for estetics. 

// experimental code that should not directly mix with production level code. 

#define STEPS_FOR_LNL 1000
#define INTEGRATION_GENERATIONS 100000
#define NR_STEPS 30


#include <sstream>
#include "priors/ExponentialPrior.hpp"
#include "proposals/BranchIntegrator.hpp"
#include "ProposalRegistry.hpp"
#include "parameters/BranchLengthsParameter.hpp"

std::pair<double,double> getMeanAndVar (const std::vector<double> &data )
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



void SampleMaster::branchLengthsIntegration()  
{
  assert(runs.size() == 1 );   
  auto &run = runs[0];   
  auto &chains = run.getChains(); 
  assert(chains.size() == 1); 
  auto &chain = chains[0]; 
  
  stringstream ss; 
  ss << cl.getRunid() << ".tree.tre" ; 

  TreePrinter tp(true, true, false);   
  TreePrinter tp2(true, true, true); 
  auto tralnPtr = chain.getTralnPtr(); 
  auto& traln  = *tralnPtr  ; 

  auto eval = chain.getEvaluator()->clone(); 

  vector<unique_ptr<AbstractParameter> > params; 
  params.emplace_back(unique_ptr<AbstractParameter>(new BranchLengthsParameter( 0 ,0))); 
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    params[0]->addPartition(i);

  double lambda   =  10 ; 

  params[0]->setPrior(make_shared< ExponentialPrior>(lambda));

  auto p = unique_ptr<BranchIntegrator>(new BranchIntegrator (ProposalRegistry::initBranchLengthMultiplier)); 
  vector<unique_ptr<AbstractProposal> >  proposals;   
  proposals.push_back( std::move(p) ); 
  proposals[0]->addPrimaryParameter( std::move(params[0])); 
  std::vector<AbstractParameter*> paramView; 
  for(auto &v : params)
    paramView.push_back(v.get());

  std::vector<ProposalSet> pSets; 

  Chain integrationChain(masterRand.generateSeed(), tralnPtr, proposals, pSets, eval->clone() );   

  auto branches =  traln.extractBranches(paramView);
  auto ps = integrationChain.getProposalView(); 
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0]); 

  integrationChain.suspend(); 

  for(auto &branch : branches)
    {      
      // setup 
      integrationChain.resume(false, false ); 
      traln.setBranch(branch, paramView); 
      eval->evaluate(traln, branch, true); 
      integrationChain.reinitPrior();

      double minHere = 1000; 
      double maxHere = 0; 
      Branch initBranch = branch; 

      traln.setBranch(branch, paramView); 
      eval->evaluate(traln, branch, true); 
      integrator->setToPropose(branch); 

      ////////////////////////
      // integrate branch   //
      ////////////////////////      
      stringstream ss; 
      ss << "samples." << cl.getRunid()<< "." << branch.getPrimNode() << "-" << branch.getSecNode()   <<  ".tab" ;
      ofstream thisOut (ss.str());       
      // run the chain to integrate 
      tout << "integrating branch " << branch << endl; 
      for(int i = 0; i < INTEGRATION_GENERATIONS; ++i) 
	{	  
	  integrationChain.step();
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
      
      // get branch lengths 
      ss.str(std::string()); 
      ss << "lnl." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str());

      tout << "evaluating branch lengths for " << branch << endl; 
      
      if(maxHere != minHere)
	{
	  for(double i = minHere; i < maxHere+0.00000001 ; i+= (maxHere-minHere)/ STEPS_FOR_LNL)
	    {
	      assert(0);
	      // double tmp = branch.getInternalLength(traln,i); 
	      Branch b = branch; 
	      // b.setLength(tmp); 
	      traln.setBranch(b, paramView[0]); 
	      
	      // traln.setBranchLengthBounded(tmp, 0, branch.findNodePtr(traln)); 
	      eval->evaluate(traln, branch, false);
	      double lnl = traln.getTr()->likelihood; 
	      
	      thisOut << i << "\t" << setprecision(std::numeric_limits<double>::digits10) << lnl << endl; 
	    }
	}
      else
	thisOut << minHere << "\t" << "NA" << endl; 
      thisOut.close(); 

      Branch tmpBranch = branch; 
      tout << "optimizing the branch using nr" << endl; 
      ss.str(std::string());
      ss << "nr-length." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str()); 

      double result = 0;  
      // double curVal = branch.getInternalLength(traln,0.1);
      double curVal = 0; 
      // TODO 
      double secDerivative = 0; 
      double firstDerivative = 0; 

      double prevVal = curVal; 
      for(int i = 0; i < NR_STEPS; ++i )
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

      // double something = tmpBranch.getInternalLength(traln, prevVal); 
      // TODO 
      double something = 0; 
      assert(0);


#if HAVE_PLL != 0
      makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
#else 
      makenewzGeneric(traln.getTr(), 
		      branch.findNodePtr(traln), branch.getInverted().findNodePtr(traln),
		      &something, 1, &result,  &firstDerivative, &secDerivative, lambda, FALSE); 
#endif
      
      thisOut << prevVal << "\t" << firstDerivative << "\t" << secDerivative << endl; 

      thisOut.close(); 
      
      // reset 
      traln.setBranch(initBranch, paramView); 
    }

  tout << "finished!" << endl; 
}

