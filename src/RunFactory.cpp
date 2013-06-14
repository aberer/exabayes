#include "RunFactory.hpp"
#include "ProposalRegistry.hpp"
#include "output.h"
#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"

// not to be confused with a fun factory...

void RunFactory::addStandardParameters(vector<RandomVariable> &vars, const TreeAln &traln )
{
  vector<bool> categoryIsActive( NUM_PROP_CATS, false );
  for(auto v : vars )
    categoryIsActive[v.getCategory()] = true; 

  nat highestId = vars.size() == 0 ? 0 : vars[vars.size()-1].getId(); 

  // add standard stuff, if not defined yet
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    {
      category_t catIter = category_t(i) ; 
      switch(catIter)
	{	  
	case TOPOLOGY: 
	  {
	    // deafult: everything linked 
	    if(not categoryIsActive[catIter])
	      {
		RandomVariable r(catIter,highestId);
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  r.addPartition(j); 
		vars.push_back(r); 
	      }
	    break; 
	  }

	case BRANCH_LENGTHS: 
	  {
	    // default: everything is linked
	    if(not categoryIsActive[catIter])
	      {
		RandomVariable r(catIter,highestId);
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  r.addPartition(j); 
		vars.push_back(r); 
	      }
	    break; 
	  }

	case FREQUENCIES: 
	  {
	    if(not categoryIsActive[catIter] )
	      {		
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  {		    
		    RandomVariable r(catIter, highestId);
		    r.addPartition(j); 
		    vars.push_back(r);	   
		  }
	      }
	    break; 
	  }

	case SUBSTITUTION_RATES: 
	  {
	    if(not categoryIsActive[catIter] )
	      {
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  {
		    RandomVariable r(catIter, highestId);
		    r.addPartition(j); 
		    vars.push_back(r);	   
		  }
	      }
	    break;  
	  }

	case RATE_HETEROGENEITY: 
	  {
	    if(not categoryIsActive[catIter])
	      {
		for(int j = 0; j < traln.getNumberOfPartitions();++ j)
		  {
		    RandomVariable r(catIter, highestId);
		    r.addPartition(j);
		    vars.push_back(r);
		  }
	      }
	    break; 
	  }
	  
	case AA_MODEL: 
	  {
	    if( not categoryIsActive[catIter])
	      {
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  {
		    pInfo* partition = traln.getPartition(j);
		    if(partition->dataType == AA_DATA)
		      {
			RandomVariable r(catIter, highestId);
			r.addPartition(j); 
			vars.push_back(r);
		      }
		  }
	      }
	    break; 
	  }
	default: assert(0); 
	}
    }
}


void RunFactory::addStandardPrior(RandomVariable &var, const TreeAln& traln )
{
  switch(var.getCategory())			// TODO such switches should be part of an object
    {
    case TOPOLOGY:  
      var.setPrior(make_shared<UniformPrior>(0,0)); // TODO : proper topology prior? 
      break; 
    case BRANCH_LENGTHS: 
      var.setPrior(shared_ptr<AbstractPrior>(new ExponentialPrior(10.0)));
      break; 
    case FREQUENCIES: 
      {
	pInfo *partition = traln.getPartition(var.getPartitions()[0]);
	assert(partition->dataType == DNA_DATA || partition->dataType == AA_DATA); 

	vector<double>badHardcoded; 
	for(int i = 0; i < partition->states; ++i)
	  badHardcoded.push_back(1.); 
	var.setPrior(shared_ptr<AbstractPrior>(new DirichletPrior(badHardcoded ))); 
      }
      break; 
    case SUBSTITUTION_RATES: 
      {
	pInfo *partition = traln.getPartition(var.getPartitions()[0]);
	assert(partition->dataType == DNA_DATA); 
	
	vector<double> subAlpha = {1,1,1,1,1,1}; 
	var.setPrior(shared_ptr<AbstractPrior>(new DirichletPrior( subAlpha ))); 
      }
      break; 
    case RATE_HETEROGENEITY: 
      var.setPrior(shared_ptr<AbstractPrior>(new UniformPrior(1e-6, 200)));     
      break; 
    case AA_MODEL : 
      assert(NOT_IMPLEMENTED); 
      break; 
    default: assert(0); 
    }
}


void RunFactory::addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<RandomVariable> &variables)
{
  auto generalPriors = priorInfo.getGeneralPriors();
  auto specificPriors = priorInfo.getSpecificPriors();

  for(RandomVariable &v : variables)
    {
      auto partitionIds = v.getPartitions(); 

      // try adding a partition specific prior 
      auto idMap = specificPriors[ v.getCategory() ]; 
      shared_ptr<AbstractPrior> thePrior = nullptr; 
      for(nat partId : partitionIds)	
	{	  	  
	  if(idMap.find(partId) != idMap.end()) // found 
	    {
	      if(thePrior != nullptr)
		{
		  cerr << "while setting prior for random variable " << v << ": it seems you have defined a prior for more than one partition of a set of partitions that is linked for a specific variable to estimate. Please only specify exactly one prior per variable for a set of linked partitions."  << endl; 
		  exit(1); 
		}
	      else 
		thePrior = idMap.at(partId); 
	    }	 	    
	}

      // use a general prior, if we have not found anything
      if(thePrior == nullptr)	  
	{
	  thePrior = generalPriors.at(v.getCategory()); 	  
	  if(thePrior != nullptr)
	    cout << "using GENERAL prior for variable " <<  endl; 
	}
      else 
	{
	  cout << "using PARTITION-SPECIFIC prior for variable " <<  endl; 
	}
      
      if(thePrior != nullptr)
	v.setPrior(thePrior);
      else 
	{	  
	  addStandardPrior(v, traln);
	  cout << "using STANDARD prior for variable "  << v << endl; 
	}
    }
}

void RunFactory::configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, vector<unique_ptr<AbstractProposal> > &proposals)
{
  randomVariables = partitionParams.getParameters();
  addStandardParameters(randomVariables, traln);
  addPriorsToVariables(traln, priorInfo, randomVariables);

  ProposalRegistry reg; 

  vector<RandomVariable> blRandVars; 
  for(auto v : randomVariables)
    if(v.getCategory() == BRANCH_LENGTHS)
      blRandVars.push_back(v); 

  for(auto v : randomVariables)
    {
      if(typeid(v.getPrior()) == typeid(shared_ptr<FixedPrior>))
	continue; 

      vector<unique_ptr<AbstractProposal> > tmpResult;  

      reg.getProposals(v.getCategory(), propConfig, tmpResult); 
      for(auto  &p : tmpResult )
	{
	  p->addPrimVar(v);
	  if(v.getCategory() == TOPOLOGY)
	    for(auto blRandVar : blRandVars)
	      p->addSecVar(blRandVar); 
	}

      // dammit...
      for(auto it = tmpResult.begin(); it != tmpResult.end(); )
	{
	  proposals.emplace_back(std::move(*it));
	  it = tmpResult.erase(it);
	}
    }


  // lets see if it worked 
  // seems to work....
  tout << endl << "RandomVariables to be integrated: " << endl; 
  for(auto v : randomVariables)
    tout << v  << endl; 
  tout << endl; 
  

  double sum = 0; 
  for(auto &p : proposals)
    sum += p->getRelativeWeight(); 
  
  tout << "Will employ the following proposal mixture (frequency,type,affected variables): " << endl; 
  for(auto &p : proposals )
    {
      tout << setprecision(2) << p->getRelativeWeight() / sum * 100 <<   "%\t" ; 
      p->printShort(tout ) ; 
      tout << endl; 
    }
  tout << setprecision(2) << fixed << endl ; 
}

