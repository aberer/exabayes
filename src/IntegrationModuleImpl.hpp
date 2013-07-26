// just spliced this off for estetics. 

// experimental code that should not directly mix with production level code. 

#define STEPS_FOR_LNL  1000 
#define INTEGRATION_GENERATIONS 10000

// #define STEPS_FOR_LNL  10 
// #define INTEGRATION_GENERATIONS 
#define NR_STEPS 30


#include "BoundsChecker.hpp"
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

  ofstream tFile( ss.str());   
  TreePrinter tp(true, true, false);   
  TreePrinter tp2(true, true, true); 
  auto tralnPtr = chain.getTralnPtr(); 
  auto& traln  = *tralnPtr  ; 
  tFile << tp.printTree(traln) << endl; 
  tFile << tp2.printTree(traln) << endl; 
  tFile.close(); 

  auto eval = chain.getEvaluator()->clone(); 

  vector<unique_ptr<AbstractParameter> > vars; 
  vars.emplace_back(unique_ptr<AbstractParameter>(new BranchLengthsParameter( 0 ))); 
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    vars[0]->addPartition(i);


  double lambda   =  10 ; 

  vars[0]->setPrior(make_shared< ExponentialPrior>(lambda));

  auto p = unique_ptr<BranchIntegrator>(new BranchIntegrator (ProposalRegistry::initBranchLengthMultiplier)); 
  vector<unique_ptr<AbstractProposal> >  proposals;   
  proposals.push_back( std::move(p) ); 
  proposals[0]->addPrimVar( std::move(vars[0])); 

  Chain integrationChain(masterRand.generateSeed(), tralnPtr, proposals, eval->clone() );   

  auto branches =  traln.extractBranches();
  auto ps = integrationChain.getProposalView(); 
  assert(ps.size() == 1 );   
  auto integrator = dynamic_cast<BranchIntegrator*>(ps[0]); 

  integrationChain.suspend(false); 

  for(auto &branch : branches)
    {      
      // setup 
      integrationChain.resume(false, false ); 
      traln.setBranch(branch); 
      eval->evaluate(traln, branch, true); 
      integrationChain.reinitPrior();

      double minHere = 1000; 
      double maxHere = 0; 
      Branch initBranch = branch; 

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
	  auto elem = traln.getBranch(branch.findNodePtr(traln)); 
	  auto iLen = elem.getInterpretedLength(traln);
	  if (i % 10 == 0)
	    thisOut << iLen << endl; 
	  if(iLen < minHere)
	    minHere = iLen; 
	  if(maxHere < iLen)
	    maxHere = iLen; 
	} 
      thisOut.close();


      /////////////////////////////////////
      // get branch length likelihoods   //
      /////////////////////////////////////
      ss.str(std::string()); 
      ss << "lnl." << cl.getRunid() << "." << branch.getPrimNode() << "-" << branch.getSecNode() << ".tab"; 
      thisOut.open(ss.str());

      tout << "evaluating branch lengths for " << branch << endl; 
      
      if(maxHere != minHere)
	{
	  for(double i = minHere; i < maxHere+0.00000001 ; i+= (maxHere-minHere)/ STEPS_FOR_LNL)
	    {
	      double tmp = branch.getInternalLength(traln,i); 
	      Branch b = branch; 
	      b.setLength(tmp); 
	      traln.setBranch(b); 

	      eval->evaluate(traln, branch, false);
	      double lnl = traln.getTr()->likelihood; 
	      
	      thisOut << i << "\t" << setprecision(std::numeric_limits<double>::digits10) << lnl << endl; 
	    }
	}
      else
	thisOut << minHere << "\t" << "NA" << endl; 
      thisOut.close(); 


      /////////////////////////
      // optimize the branch //
      /////////////////////////
      Branch tmpBranch = branch; 
      tout << "optimizing the branch using nr" << endl; 
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


      ///////////////////////////
      // BRANCH CHANGE IMPACT  //
      ///////////////////////////
      ss.str(""); 
      ss << "branchImpact." << cl.getRunid()<< "." << tmpBranch.getPrimNode() << "-" << tmpBranch.getSecNode()   <<  ".tab" ;
      thisOut.open(ss.str()); 
      
      tmpBranch.setLength(result);       
      traln.setBranch(tmpBranch);
      
      std::vector<Branch> firstOrderDesc; 
      if(not traln.isTipNode(tmpBranch.findNodePtr(traln) ))
	{
	  auto res = traln.getDescendants(tmpBranch); 
	  firstOrderDesc.push_back(res.first.getInverted()); 
	  firstOrderDesc.push_back(res.second.getInverted()); 
	}
      if(not traln.isTipNode(tmpBranch.getInverted().findNodePtr(traln)) )
	{
	  auto res = traln.getDescendants(tmpBranch.getInverted()); 
	  firstOrderDesc.push_back(res.first.getInverted()); 	  
	  firstOrderDesc.push_back(res.second.getInverted()); 	  
	}

      std::vector<Branch> secOrderDesc; 
      for(auto &b : firstOrderDesc)
	{
	  if(not b.isTipBranch(traln))
	    {
	      auto res = traln.getDescendants(b);
	      secOrderDesc.push_back(res.first.getInverted()); 
	      secOrderDesc.push_back(res.second.getInverted()); 
	    }
	}

      // get the real branch lengths 
      for(auto &b : firstOrderDesc)
	b = traln.getBranch(b); 
      for(auto &b : secOrderDesc)
	b = traln.getBranch(b);       

      std::vector<Branch> allDescendants; 
      allDescendants.insert(allDescendants.end(), firstOrderDesc.begin(), firstOrderDesc.end()); 
      allDescendants.insert(allDescendants.end(), secOrderDesc.begin(), secOrderDesc.end()); 

      // meaning of vector: 
      // (mean_opt, var_opt), (mean_half, var_half), (mean_double, var_double)
      std::unordered_map<Branch,std::vector<double>,BranchHashNoLength, BranchEqualNoLength> branch2sampleStat; 
      for(auto &b : allDescendants)
	branch2sampleStat[b]= std::vector<double>(); 
      
      // try out these branch lengths for the main branch 
      for(auto alt: {tmpBranch.getLength(), sqrt(tmpBranch.getLength()), pow(tmpBranch.getLength(),2) })
	{      
	  // tout << "tree is " << traln << std::endl; 
	  tmpBranch.setLength(alt) ; 
	  if(not BoundsChecker::checkBranch(tmpBranch))
	    BoundsChecker::correctBranch(tmpBranch); 
	  // tout << "setting branch under examination to " << tmpBranch << std::endl; 
	  
	  traln.setBranch(tmpBranch); 
	  
	  for(const Branch &b :allDescendants)
	    {
	      traln.setBranch(b); 
	
	      eval->evaluate(traln,b, true); 
	      integrationChain.reinitPrior(); 

	      std::vector<double> samples; 

	      integrator->setToPropose(b);
	      for(int i = 0; i < INTEGRATION_GENERATIONS; ++i )
		{
		  integrationChain.step();
		  auto result = traln.getBranch(b); 
		  auto iLen = result.getInterpretedLength(traln);
		  if(i  % 10  == 0)
		    samples.push_back(iLen);
		}
	  
	      auto res = getMeanAndVar(samples); 	  
	      branch2sampleStat[b].push_back(res.first); 
	      branch2sampleStat[b].push_back(res.second); 
	  
	      traln.setBranch(b);
	    }
	}

      Branch::printLength = false;  
      for(auto &f : firstOrderDesc )
	{
	  assert(branch2sampleStat[f].size() == 6 ); 
	  assert(branch2sampleStat.find(f) != branch2sampleStat.end());
	  thisOut << "FIRST\t" << f << "\tOPT\t"  <<  branch2sampleStat[f][0]  << "\t" << branch2sampleStat[f][1] << std::endl; 
	  thisOut << "FIRST\t" << f << "\tHALF\t"  << branch2sampleStat[f][2]  << "\t" << branch2sampleStat[f][3] << std::endl; 
	  thisOut << "FIRST\t" << f << "\tDOUBLE\t"  << branch2sampleStat[f][4]  << "\t" << branch2sampleStat[f][5] << std::endl; 
	}
      for(auto &s : secOrderDesc )
	{
	  assert(branch2sampleStat[s].size() == 6 ); 
	  assert(branch2sampleStat.find(s) != branch2sampleStat.end());
	  thisOut << "SECOND\t" << s << "\tOPT\t"  << branch2sampleStat[s][0]  << "\t" << branch2sampleStat[s][1] << std::endl; 
	  thisOut << "SECOND\t" << s << "\tHALF\t"  << branch2sampleStat[s][2]  << "\t" << branch2sampleStat[s][3] << std::endl; 
	  thisOut << "SECOND\t" << s << "\tDOUBLE\t"  << branch2sampleStat[s][4]  << "\t" << branch2sampleStat[s][5] << std::endl; 
	}
      thisOut.close();
      
      Branch::printLength = true; 
      
      // reset 
      traln.setBranch(initBranch); 
      tout << std::endl; 
    }

  tout << "finished!" << endl; 
}

