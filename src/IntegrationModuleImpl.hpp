// just spliced this off for estetics. 

// experimental code that should not directly mix with production level code. 

#define STEPS_FOR_LNL 1000
#define INTEGRATION_GENERATIONS 10000
#define NR_STEPS 30

#include <sstream>
#include "proposals/BranchIntegrator.hpp"
#include "ProposalRegistry.hpp"

void SampleMaster::branchLengthsIntegration()  
{
  assert(runs.size() == 1 );   
  auto &run = runs[0];   
  auto &chains = run.getChains(); 
  assert(chains.size() == 1); 
  auto &chain = chains[0]; 
  
  stringstream ss; 
  ss << cl.getRunid() << ".tree.tre" ; 

  ofstream tFile( ss.str());   
  TreePrinter tp(true, true, false);   
  TreePrinter tp2(true, true, true); 
  auto tralnPtr = chain.getTralnPtr(); 
  auto& traln  = *tralnPtr  ; 
  tFile << tp.printTree(traln) << endl; 
  tFile << tp2.printTree(traln) << endl; 
  tFile.close(); 
  
  auto eval = chain.getEvaluatorPtr();

  auto r = make_shared<RandomVariable>(Category::BRANCH_LENGTHS, 0);
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    r->addPartition(i);

  vector<shared_ptr<RandomVariable> > vars =  {r} ; 

  double lambda   =  10 ; 

  vars[0]->setPrior(make_shared< ExponentialPrior>(lambda));

  auto p = unique_ptr<BranchIntegrator>(new BranchIntegrator (ProposalRegistry::initBranchLengthMultiplier)); 
  vector<unique_ptr<AbstractProposal> >  proposals;   
  proposals.push_back( std::move(p) ); 
  proposals[0]->addPrimVar(vars[0]); 

  Chain integrationChain(masterRand.generateSeed(), tralnPtr, proposals, eval );   

  auto branches =  traln.extractBranches();
  auto ps = integrationChain.getProposalView(); 
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0]); 

  for(auto &branch : branches)
    {      
      double minHere = 1000; 
      double maxHere = 0; 

      traln.setBranch(branch); 
      integrator->setToPropose(branch); 
      
      stringstream ss; 
      ss << "samples." << cl.getRunid()<< "." << branch.getPrimNode() << "-" << branch.getSecNode()   <<  ".tab" ;
      ofstream thisOut (ss.str()); 
      
      // run the chain to integrate 
      cout << "integrating branch " << branch << endl; 
      for(int i = 0; i < INTEGRATION_GENERATIONS; ++i) 
	{	  
	  integrationChain.step();
	  traln.setBranch(branch); 	  

	  double length = branch.getInterpretedLength(traln); 
	  thisOut << length << endl; 
	  
	  if(length < minHere)
	    minHere = length; 
	  if(maxHere < length )
	    maxHere = length; 
	}      

      thisOut.close();
      
      // get branch lengths 
      ss.str(std::string()); 
      ss << "lnl." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str());

      cout << "evaluating branch lengths for " << branch << endl; 
      
      if(maxHere != minHere)
	{
	  for(double i = minHere; i < maxHere+0.00000001 ; i+= (maxHere-minHere)/ STEPS_FOR_LNL)
	    {
	      double tmp = branch.getInternalLength(traln,i); 
	      Branch b = branch; 
	      b.setLength(tmp); 
	      traln.setBranch(b); 
	      
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
      cout << "optimizing the branch using nr" << endl; 
      ss.str(std::string());
      ss << "nr-length." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str()); 

      double result = 0;  
      double curVal = branch.getInternalLength(traln,0.1); 
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
	  tmpBranch.setLength(result);
	  thisOut << prevVal <<  "\t" << firstDerivative << "\t" << secDerivative << endl; 	
	  prevVal = tmpBranch.getInterpretedLength(traln); 
	  curVal = result; 
	} 

      double something = tmpBranch.getInternalLength(traln, prevVal); 

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
      traln.setBranch(branch); 
    }

  cout << "finished!" << endl; 
}
