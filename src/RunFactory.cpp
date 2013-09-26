#include <set>
#include <iostream>
#include <unordered_map>

#include "RunFactory.hpp"
#include "ProposalRegistry.hpp"

#include "proposers/AbstractProposer.hpp"
#include "proposers/MultiplierProposal.hpp"
#include "proposers/DirichletProposal.hpp"
#include "proposers/SlidingProposal.hpp"
#include "ParallelSetup.hpp"

#include "parameters/TopologyParameter.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/FrequencyParameter.hpp"
#include "parameters/RevMatParameter.hpp"
#include "parameters/RateHetParameter.hpp"

#include "priors/AbstractPrior.hpp"
#include "priors/ExponentialPrior.hpp"
#include "priors/UniformPrior.hpp"
#include "priors/DirichletPrior.hpp"
#include "priors/FixedPrior.hpp"


void RunFactory::addStandardParameters(vector<unique_ptr<AbstractParameter> > &vars, const TreeAln &traln )
{
  auto categories =  std::set<Category>{}; 

  for(auto &v : vars)
    categories.insert(v->getCategory()); 

  int highestId = -1; 
  for(auto &v : vars )
    {
      int id = v->getId(); 
      if(highestId < id)
	highestId = id; 
    }
  ++highestId;

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
	    auto r = CategoryFuns::getParameterFromCategory(catIter, highestId, 0 );
	    ++highestId; 
	    for(nat j = 0; j < traln.getNumberOfPartitions(); ++j)
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
	    for(nat j = 0; j < traln.getNumberOfPartitions(); ++j)
	      {		
		auto r = CategoryFuns::getParameterFromCategory(catIter, highestId,j) ; 
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
      var->setPrior( std::make_shared<UniformPrior>(0,0) ); // TODO : proper topology prior? 
      break; 
    case Category::BRANCH_LENGTHS: 
      var->setPrior( std::make_shared<ExponentialPrior>(10.0) );
      break; 
    case Category::FREQUENCIES: 
      {
	pInfo *partition = traln.getPartition(var->getPartitions()[0]);
	assert(partition->dataType == DNA_DATA || partition->dataType == AA_DATA); 
	var->setPrior(make_shared<DirichletPrior>(vector<double>(partition->states , 1.))); 
      }
      break; 
    case Category::SUBSTITUTION_RATES: 
      {
	pInfo *partition = traln.getPartition(var->getPartitions()[0]);	;
	var->setPrior(make_shared<DirichletPrior>( vector<double> (numStateToNumInTriangleMatrix(partition->states), 1.) )); 
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
		  ParallelSetup::genericExit(-1); 
		}
	      else 
		thePrior = idMap.at(partId); 
	    }	 	    
	}

      // use a general prior, if we have not found anything
      if(thePrior == nullptr)	  
	{
	  thePrior = generalPriors.at(int(v->getCategory())); // BAD
	  // if(thePrior != nullptr)
	  //   cout << "using GENERAL prior for variable " <<  endl; 
	}
      else 
	{
	  // tout << "using PARTITION-SPECIFIC prior for variable " <<  endl; 
	}
      
      if(thePrior != nullptr)
	v->setPrior(thePrior);
      else 
	{	  
	  addStandardPrior(v.get(), traln);
	  // tout << "using STANDARD prior for variable "  << v.get() << endl; 
	}
    }
}




void RunFactory::addSecondaryParameters(AbstractProposal* proposal, const std::vector<unique_ptr<AbstractParameter> > &allParameters)
{
  // get all branch length pararemeters   
  std::vector<unique_ptr<AbstractParameter> > blParameters; 
  for(auto &v : allParameters)
    if(v->getCategory() == Category::BRANCH_LENGTHS)
      blParameters.push_back(std::unique_ptr<AbstractParameter>(v->clone())); 

  bool needsBl = false; 
  std::unordered_set<nat> myPartitions; 
  for(auto &v: proposal->getPrimaryParameterView())
    {
      needsBl |= ( v->getCategory() == Category::FREQUENCIES 
		   || v->getCategory() == Category::SUBSTITUTION_RATES
		   || v->getCategory() == Category::TOPOLOGY
		   ); 
      auto ps = v->getPartitions();
      myPartitions.insert(ps.begin(), ps.end()); 
    }

  if(needsBl)
    {
      for(auto &blRandVar : blParameters)
	{
	  auto partitions =  blRandVar->getPartitions();
	  if(std::any_of(partitions.begin(),partitions.end(), 
			 [&](nat i){ return myPartitions.find(i) != myPartitions.end(); }))
	    proposal->addSecondaryParameter(std::unique_ptr<AbstractParameter>(blRandVar->clone())); 
	}
    }
}


auto RunFactory::produceProposals(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, bool componentWiseMH, std::vector<ProposalSet> &resultPropSet)
  -> std::vector<std::unique_ptr<AbstractProposal> > 
{
  std::vector<std::unique_ptr<AbstractProposal> > proposals; 

  randomVariables = partitionParams.getParameters(); 
  addStandardParameters(randomVariables, traln);  
  addPriorsToVariables(traln, priorInfo, randomVariables);

  ProposalRegistry reg; 

  // instantiate all proposals that integrate over one parameter 
  for(auto &v : randomVariables)
    {
      if(dynamic_cast<FixedPrior*>(v->getPrior()) != nullptr ) // is it a fixed prior?
	continue;
      
      // TODO  aa-models as well later 
      if( componentWiseMH && 
	 ( v->getCategory() == Category::FREQUENCIES
	   || v->getCategory() == Category::BRANCH_LENGTHS
	   || v->getCategory() == Category::SUBSTITUTION_RATES
	   || v->getCategory() == Category::RATE_HETEROGENEITY))
	continue; 
	
      std::vector<std::unique_ptr<AbstractProposal> > 
	tmpResult = reg.getSingleParameterProposals(v->getCategory(), propConfig,  traln); 
      for(auto &p : tmpResult )
	{
	  p->addPrimaryParameter(std::unique_ptr<AbstractParameter>(v->clone()));
	  addSecondaryParameters(p.get(), randomVariables); 
	}

      // dammit...
      for(auto it = tmpResult.begin(); it != tmpResult.end(); )
	{
	  proposals.emplace_back(std::move(*it));
	  it = tmpResult.erase(it);
	}
    }

  // instantiate proposals that integrate over multiple over an entire
  // category gather all parameters that we can integrate over
  // together in a partitioned manner and output a good set of
  // proposals
  nat blCtr = 0; 
  std::vector<AbstractParameter*> mashableParameters; 
  for(auto &v : randomVariables)
    {
      if(dynamic_cast<FixedPrior*>(v->getPrior()) != nullptr ) // is it a fixed prior?
	continue;

      // those are prototypes! non-owning pointers 
      switch(v->getCategory())
	{
	case Category::SUBSTITUTION_RATES: 
	case Category::FREQUENCIES:
	case Category::RATE_HETEROGENEITY: 
	  mashableParameters.push_back(v.get()); 
	  break; 
	case Category::BRANCH_LENGTHS:
	  {
	    mashableParameters.push_back(v.get()); 
	    ++blCtr; 
	  }
	default: 
	  ;
	}
    }  
  
#ifdef UNSURE
  // what about fixed priors? 
  assert(0); 
#endif
  // TODO that's all a bit cumbersome, has to be re-designed a bit 

  if(componentWiseMH)
    {
      auto cat2param = std::unordered_map<Category, std::vector<AbstractParameter*>, CategoryHash >{}; 
      for(auto &p:  mashableParameters)
	cat2param[p->getCategory()].push_back(p); 

      for(auto &elem : cat2param)
	{
	  auto proposalsForSet = reg.getSingleParameterProposals(elem.first, propConfig, traln);
	  
	  // this and the above for-loop essentially produce all
	  // proposals, that we'd also obtain in the default case =>
	  // but we need a set of it
	  for(auto &proposalType : proposalsForSet)
	    {
	      auto lP = std::vector<std::unique_ptr<AbstractProposal> >{}; 	      
	      for(auto &p : elem.second)
		{
		  auto proposalClone = std::unique_ptr<AbstractProposal>(proposalType->clone()); 
		  proposalClone->addPrimaryParameter(std::unique_ptr<AbstractParameter>(p->clone()));
		  addSecondaryParameters(proposalClone.get(), randomVariables); 
		  lP.push_back(std::move(proposalClone));		  
		} 
	      proposalType->setInSetExecution(true);

	      resultPropSet.emplace_back(lP[0]->getRelativeWeight(), std::move(lP)); 
	    }
	}
    }

  return proposals; 
}


std::vector<std::unique_ptr<AbstractParameter> > RunFactory::getRandomVariables() const 
{
  vector<unique_ptr<AbstractParameter> > result; 
  for(auto &r : randomVariables)
    result.push_back(std::unique_ptr<AbstractParameter>(r->clone()));
  return result; 
}
