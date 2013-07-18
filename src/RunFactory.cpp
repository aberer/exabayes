#include <set>

#include "RunFactory.hpp"
#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"

#include "parameters/TopologyParameter.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/FrequencyParameter.hpp"
#include "parameters/RevMatParameter.hpp"
#include "parameters/RateHetParameter.hpp"



// not to be confused with a fun factory...
void RunFactory::addStandardParameters(vector<unique_ptr<AbstractParameter> > &vars, const TreeAln &traln )
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
	  // force to have everything linked with those 
	case Category::TOPOLOGY: 
	case Category::BRANCH_LENGTHS: 
	  {
	    auto r = CategoryFuns::getParameterFromCategory(catIter, highestId );
	    ++highestId; 
	    for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
	      r->addPartition(j); 
	    vars.push_back(std::move(r)); 
	  }
	  break; 
	case Category::AA_MODEL:	
	  std::cout << "TODO define,  when AA_MODLE should be added <= runfactory.cpp" << std::endl; 
	  break; 	  
	  // a new parameter per partition 
	case Category::SUBSTITUTION_RATES: 
	case Category::RATE_HETEROGENEITY:
	case Category::FREQUENCIES:
	  {
	    for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
	      {
		auto r = CategoryFuns::getParameterFromCategory(catIter, highestId) ; 
		++highestId; 
		r->addPartition(j); 
		vars.push_back(std::move(r)); 
	      }	    
	  }
	  break; 
	default: assert(0); 
	}
    }
}


void RunFactory::addStandardPrior(AbstractParameter* var, const TreeAln& traln )
{
  switch(var->getCategory())			// TODO such switches should be part of an object
    {
    case Category::TOPOLOGY:  
      var->setPrior( make_shared<UniformPrior>(0,0)); // TODO : proper topology prior? 
      break; 
    case Category::BRANCH_LENGTHS: 
      var->setPrior(make_shared<ExponentialPrior>(10.0));
      break; 
    case Category::FREQUENCIES: 
      {
	pInfo *partition = traln.getPartition(var->getPartitions()[0]);
	assert(partition->dataType == DNA_DATA || partition->dataType == AA_DATA); 

	vector<double>badHardcoded; 
	for(int i = 0; i < partition->states; ++i)
	  badHardcoded.push_back(1.); 
	var->setPrior(make_shared<DirichletPrior>(badHardcoded)); 
      }
      break; 
    case Category::SUBSTITUTION_RATES: 
      {
	pInfo *partition = traln.getPartition(var->getPartitions()[0]);
	assert(partition->dataType == DNA_DATA); 
	
	vector<double> subAlpha = {1,1,1,1,1,1}; 
	var->setPrior(make_shared<DirichletPrior>( subAlpha )); 
      }
      break; 
    case Category::RATE_HETEROGENEITY: 
      var->setPrior(make_shared<UniformPrior>(1e-6, 200));     
      break; 
    case Category::AA_MODEL : 
      assert(NOT_IMPLEMENTED); 
      break; 
    default: assert(0); 
    }
}


void RunFactory::addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<unique_ptr<AbstractParameter> > &variables)
{
  auto generalPriors = priorInfo.getGeneralPriors();
  auto specificPriors = priorInfo.getSpecificPriors();

  for(auto &v : variables)
    {
      auto partitionIds = v->getPartitions(); 

      // try adding a partition specific prior 
      auto idMap = specificPriors[ int(v->getCategory()) ]; 
      shared_ptr<AbstractPrior> thePrior = nullptr; 
      for(nat partId : partitionIds)	
	{	  	  
	  if(idMap.find(partId) != idMap.end()) // found 
	    {
	      if(thePrior != nullptr)
		{
		  cerr << "while setting prior for random variable " << v.get() << ": it seems you have defined a prior for more than one partition of a set of partitions that is linked for a specific variable to estimate. Please only specify exactly one prior per variable for a set of linked partitions."  << endl; 
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
	  tout << "using PARTITION-SPECIFIC prior for variable " <<  endl; 
	}
      
      if(thePrior != nullptr)
	v->setPrior(thePrior);
      else 
	{	  
	  addStandardPrior(v.get(), traln);
	  tout << "using STANDARD prior for variable "  << v.get() << endl; 
	}
    }
}

void RunFactory::configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, vector<unique_ptr<AbstractProposal> > &proposals, shared_ptr<LikelihoodEvaluator> eval )
{
  addStandardParameters(randomVariables, traln);  
  addPriorsToVariables(traln, priorInfo, randomVariables);

  ProposalRegistry reg; 

  std::vector<unique_ptr<AbstractParameter> > blRandVars; 
  for(auto &v : randomVariables)
    if(v->getCategory() == Category::BRANCH_LENGTHS)
      blRandVars.push_back(std::unique_ptr<AbstractParameter>(v->clone())); 

  for(auto &v : randomVariables)
    {
      if(typeid(v->getPrior()) == typeid(shared_ptr<FixedPrior>))
	continue;

      vector<unique_ptr<AbstractProposal> > tmpResult;  

      reg.getProposals(v->getCategory(), propConfig, tmpResult, traln, eval); 
      for(auto  &p : tmpResult )
	{
	  p->addPrimVar(std::unique_ptr<AbstractParameter>(v->clone()));
	  if(v->getCategory() == Category::TOPOLOGY)
	    {
	      for(auto &blRandVar : blRandVars)
		p->addSecVar(std::unique_ptr<AbstractParameter>(blRandVar->clone())); 
	    }
	  else if(v->getCategory() == Category::FREQUENCIES || v->getCategory() == Category::SUBSTITUTION_RATES)
	    {
	      for(auto &blRandVar : blRandVars)
		p->addSecVar(std::unique_ptr<AbstractParameter>(blRandVar->clone())) ; 
	    }	  
	}

      // dammit...
      for(auto it = tmpResult.begin(); it != tmpResult.end(); )
	{
	  proposals.emplace_back(std::move(*it));
	  it = tmpResult.erase(it);
	}
    }

  tout << endl << "Parameters to be integrated: " << endl; 
  for(auto &v : randomVariables)
    tout << v.get()  << endl; 
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

