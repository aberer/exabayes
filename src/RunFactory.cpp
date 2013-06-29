#include <set>

#include "RunFactory.hpp"
#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"



// not to be confused with a fun factory...

void RunFactory::addStandardParameters(vector<RandomVariablePtr> &vars, const TreeAln &traln )
{
  std::set<Category> categories; 

  for(auto &v : vars)
    categories.insert(v->getCategory()); 

  nat highestId = vars.size() == 0 ? 0 : vars[vars.size()-1]->getId(); 

  // add standard stuff, if not defined yet
  for(auto &cat : CategoryFuns::getAllCategories())
    {
      Category  catIter = cat; 
      if(categories.find(cat) != categories.end())
	continue; 

      switch(catIter)
	{	  
	case Category::TOPOLOGY: 
	  {
	    // deafult: everything linked 
	    RandomVariablePtr r(new RandomVariable(catIter,highestId));
	    for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
	      r->addPartition(j); 
	    vars.push_back(r); 
	    break; 
	  }

	case Category::BRANCH_LENGTHS: 
	  {
	    // default: everything is linked
	    RandomVariablePtr r(new RandomVariable(catIter,highestId));
	    for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
	      r->addPartition(j); 
	    vars.push_back(r); 
	    break; 
	  }

	case Category::AA_MODEL: 
	case Category::SUBSTITUTION_RATES: 
	case Category::FREQUENCIES: 
	  {
	    for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
	      {		    
		RandomVariablePtr r(new RandomVariable(catIter, highestId));
		r->addPartition(j); 
		vars.push_back(r);	   
	      }

	    break; 
	  }
	case Category::RATE_HETEROGENEITY: 
	  {
	    for(int j = 0; j < traln.getNumberOfPartitions();++ j)
	      {
		RandomVariablePtr r(new RandomVariable(catIter, highestId));
		r->addPartition(j);
		vars.push_back(r);
	      }

	    break; 
	  }
	default: assert(0); 
	}
    }
}


void RunFactory::addStandardPrior(RandomVariablePtr &var, const TreeAln& traln )
{
  switch(var->getCategory())			// TODO such switches should be part of an object
    {
    case Category::TOPOLOGY:  
      var->setPrior( PriorPtr(new UniformPrior(0,0))); // TODO : proper topology prior? 
      break; 
    case Category::BRANCH_LENGTHS: 
      var->setPrior(PriorPtr(new ExponentialPrior(10.0)));
      break; 
    case Category::FREQUENCIES: 
      {
	pInfo *partition = traln.getPartition(var->getPartitions()[0]);
	assert(partition->dataType == DNA_DATA || partition->dataType == AA_DATA); 

	vector<double>badHardcoded; 
	for(int i = 0; i < partition->states; ++i)
	  badHardcoded.push_back(1.); 
	var->setPrior(PriorPtr(new DirichletPrior(badHardcoded ))); 
      }
      break; 
    case Category::SUBSTITUTION_RATES: 
      {
	pInfo *partition = traln.getPartition(var->getPartitions()[0]);
	assert(partition->dataType == DNA_DATA); 
	
	vector<double> subAlpha = {1,1,1,1,1,1}; 
	var->setPrior(PriorPtr(new DirichletPrior( subAlpha ))); 
      }
      break; 
    case Category::RATE_HETEROGENEITY: 
      var->setPrior(PriorPtr(new UniformPrior(1e-6, 200)));     
      break; 
    case Category::AA_MODEL : 
      assert(NOT_IMPLEMENTED); 
      break; 
    default: assert(0); 
    }
}


void RunFactory::addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<RandomVariablePtr> &variables)
{
  auto generalPriors = priorInfo.getGeneralPriors();
  auto specificPriors = priorInfo.getSpecificPriors();

  for(auto &v : variables)
    {
      auto partitionIds = v->getPartitions(); 

      // try adding a partition specific prior 
      auto idMap = specificPriors[ int(v->getCategory()) ]; 
      PriorPtr thePrior = nullptr; 
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
	  thePrior = generalPriors.at(int(v->getCategory())); // BAD
	  if(thePrior != nullptr)
	    cout << "using GENERAL prior for variable " <<  endl; 
	}
      else 
	{
	  cout << "using PARTITION-SPECIFIC prior for variable " <<  endl; 
	}
      
      if(thePrior != nullptr)
	v->setPrior(thePrior);
      else 
	{	  
	  addStandardPrior(v, traln);
	  cout << "using STANDARD prior for variable "  << *v << endl; 
	}
    }
}

void RunFactory::configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, vector<ProposalPtr> &proposals, LikelihoodEvaluatorPtr &eval )
{
  randomVariables = partitionParams.getParameters();  
  addStandardParameters(randomVariables, traln);
  addPriorsToVariables(traln, priorInfo, randomVariables);

  ProposalRegistry reg; 

  vector<RandomVariablePtr> blRandVars; 
  for(auto &v : randomVariables)
    if(v->getCategory() == Category::BRANCH_LENGTHS)
      blRandVars.push_back(v); 

  for(auto v : randomVariables)
    {
      if(typeid(v->getPrior()) == typeid(shared_ptr<FixedPrior>))
	continue;

      vector<ProposalPtr> tmpResult;  

      reg.getProposals(v->getCategory(), propConfig, tmpResult, traln, eval); 
      for(auto  &p : tmpResult )
	{
	  p->addPrimVar(v);
	  if(v->getCategory() == Category::TOPOLOGY)
	    {
	      for(auto blRandVar : blRandVars)
		p->addSecVar(blRandVar); 
	    }
	  else if(v->getCategory() == Category::FREQUENCIES || v->getCategory() == Category::SUBSTITUTION_RATES)
	    {
	      for(auto blRandVar : blRandVars)
		p->addSecVar(blRandVar); 
	    }	  
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
  for(auto &v : randomVariables)
    tout << *v  << endl; 
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

