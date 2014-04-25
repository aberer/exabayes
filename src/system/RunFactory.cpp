#include <set>
#include <iostream>

#include "system/extensions.hpp"
#include "system/RunFactory.hpp"
#include "ProposalRegistry.hpp"

#include "proposers/AbstractProposer.hpp"
#include "proposers/MultiplierProposer.hpp"
#include "proposers/DirichletProposer.hpp"
#include "proposers/SlidingProposer.hpp"

#include "parameters/TopologyParameter.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/FrequencyParameter.hpp"
#include "parameters/RevMatParameter.hpp"
#include "parameters/RateHetParameter.hpp"

#include "priors/DiscreteModelPrior.hpp"
#include "priors/AbstractPrior.hpp"
#include "priors/ExponentialPrior.hpp"
#include "priors/UniformPrior.hpp"
#include "priors/DirichletPrior.hpp"

void RunFactory::addStandardParameters(std::vector<std::unique_ptr<AbstractParameter> > &params, const TreeAln &traln ) const 
{
  int highestParamId = -1; 
  if(params.size() > 0 )
    {
      auto ids = std::vector<nat>{};
      for(auto &param : params )
	ids.push_back(param->getId() ); 
      highestParamId = * (std::max_element(ids.begin(), ids.end())) ; 
    } 

  auto cat2partsUsed = std::unordered_map<Category, std::vector<bool> >{}; 
  auto cat2idOfItsKind = std::unordered_map<Category,int>{}; 
  for( auto c : CategoryFuns::getAllCategories())
    {
      cat2partsUsed[c] = std::vector<bool>(traln.getNumberOfPartitions(), false);
      cat2idOfItsKind[c] = -1; 
    }

  for(const auto &param : params)
    {
      auto cat =param->getCategory(); 
      int id = param->getIdOfMyKind(); 
      if( cat2idOfItsKind[cat] < id )
	cat2idOfItsKind[cat] = id; 
      for(auto &p : param->getPartitions())
	cat2partsUsed[cat].at(p) = true; 
    }

  // add frequency parameters to aa-revmats, if not already there 
  auto relCat = Category::FREQUENCIES; 
  auto&& paramsToAdd= std::vector<std::unique_ptr<AbstractParameter>> {};
  for( const auto &p : params)
    {
      assert(p->getPartitions().size() > 0 );
      if ( p->getCategory() == Category::SUBSTITUTION_RATES
	   &&  traln.getPartition(p->getPartitions()[0]).getDataType() == PLL_AA_DATA ) 
	{
	  for(auto part : p->getPartitions())
	    {
	      if(not cat2partsUsed.at(relCat).at(part) )
		{
		  ++highestParamId;   
		  ++cat2idOfItsKind[relCat]; 
		  paramsToAdd.push_back(CategoryFuns::getParameterFromCategory(relCat, highestParamId, cat2idOfItsKind[relCat], {part}));
		  cat2partsUsed.at(relCat).at(part) = true; 
		}
	    }
	}
    }
  
  for(auto &p : paramsToAdd)
    params.push_back(std::move(p)); 

  // add standard stuff, if not defined yet
  for(auto cat : CategoryFuns::getAllCategories())
    {
      // determine unused partitions 
      auto partsUnused = std::vector<nat>{}; 
      switch(cat)
	{	  
	  // force to have everything linked with those 
	case Category::TOPOLOGY: 
	  {
	    for(nat i = 0; i < traln.getNumberOfPartitions() ; ++ i)
	      {
		if(not cat2partsUsed.at(cat).at(i))
		  partsUnused.push_back(i);
	      }
	  }
	  break; 
	case Category::AA_MODEL:	
	  {
	    for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
	      {
		if( traln.getPartition(i).getDataType() == PLL_AA_DATA
		    && not cat2partsUsed[Category::SUBSTITUTION_RATES][i]
		    && not cat2partsUsed[Category::AA_MODEL][i] )
		  partsUnused.push_back(i);
	      }
	  } 
	  break; 
	case Category::BRANCH_LENGTHS: 
	case Category::RATE_HETEROGENEITY:
	  {
	    for(nat j = 0; j < traln.getNumberOfPartitions(); ++j)
	      {
		if(not cat2partsUsed[cat][j])
		  partsUnused.push_back(j);
	      }
	  }
	  break; 
	case Category::FREQUENCIES:
	case Category::SUBSTITUTION_RATES: 
	  {
	    for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i )
	      {
		if(not cat2partsUsed.at(cat).at(i) // not already there 
		   && traln.getPartition(i).getDataType() != PLL_AA_DATA) // we use an AA_MODEL as default for proteins 
		  partsUnused.push_back(i); 
	      }
	  }
	  break; 
	default: 
	  assert(0); 
	}

      // create parameters for unused partitions 
      switch(cat)
	{
	case Category::TOPOLOGY:
	  {
	    ++highestParamId; 
	    ++cat2idOfItsKind[cat]; 
	    params.push_back(CategoryFuns::getParameterFromCategory(cat, highestParamId, cat2idOfItsKind[cat] , partsUnused)); 
	  }
	  break; 
	case Category::FREQUENCIES: 
	case Category::AA_MODEL:	
	case Category::SUBSTITUTION_RATES: 
	case Category::BRANCH_LENGTHS: 
	case Category::RATE_HETEROGENEITY:
	  {
	    for(auto p : partsUnused)
	      {
		++highestParamId;
		++cat2idOfItsKind[cat];
		params.push_back(CategoryFuns::getParameterFromCategory(cat, highestParamId, cat2idOfItsKind[cat], {p})); 
	      }
	  }
	  break; 
	default : 
	  assert(0); 
	}
    }
}


void RunFactory::addStandardPrior(AbstractParameter* var, const TreeAln& traln )
{
  switch(var->getCategory())			// TODO such switches should be part of an object
    {
    case Category::TOPOLOGY:  
      var->setPrior( std::unique_ptr<AbstractPrior>(new UniformPrior(0,0)) ); // TODO : proper topology prior? 
      break; 
    case Category::BRANCH_LENGTHS: 
      var->setPrior( std::unique_ptr<AbstractPrior>(new ExponentialPrior(10.0) ));
      break; 
    case Category::FREQUENCIES: 
      {
	auto& partition = traln.getPartition(var->getPartitions()[0]);
	assert(partition.getDataType() == PLL_DNA_DATA || partition.getDataType() == PLL_AA_DATA); 
	var->setPrior(std::unique_ptr<AbstractPrior>(new DirichletPrior(std::vector<double>(partition.getStates() , 1.)))); 
      }
      break; 
    case Category::SUBSTITUTION_RATES: 
      {
	auto& partition = traln.getPartition(var->getPartitions()[0]);	;
	var->setPrior( make_unique<DirichletPrior>( std::vector<double>(RateHelper::numStateToNumInTriangleMatrix(partition.getStates()), 1.) )); 
      }
      break; 
    case Category::RATE_HETEROGENEITY: 
      var->setPrior(std::unique_ptr<AbstractPrior>(new UniformPrior(1e-6, 200))); 
      break; 
    case Category::AA_MODEL : 
      {
	auto modelProbs = std::unordered_map<ProtModel,double>{}; 
	for(auto model : ProtModelFun::getAllModels())
	  modelProbs[model] = 1.; 
	auto prior =  std::unique_ptr<AbstractPrior>(new DiscreteModelPrior(modelProbs));
	var->setPrior(prior);
      }
      break; 
    default: 
      assert(0); 
    }
}


void RunFactory::addPriorsToParameters(const TreeAln &traln,  const BlockPrior &priorInfo, std::vector<unique_ptr<AbstractParameter> > &variables)
{
  auto& priors = priorInfo.getPriors();

  for(auto &v : variables)
    {
      auto partitionIds = v->getPartitions(); 
      
      auto cat = v->getCategory();
      auto foundPrior = bool{false}; 

      // try to find a partition specific prior 
      for(auto iter = priors.find(cat) ; iter != priors.end() ; ++iter)
	{
	  auto &partitionsOfPrior = std::get<0>(iter->second); 
	  foundPrior = std::any_of(begin(partitionIds), end(partitionIds), [&]( nat tmp ){  return partitionsOfPrior.find(tmp) != partitionsOfPrior.end(); } ); 
	  
	  auto& prior = *(std::get<1>(iter->second).get()); 

	  if(foundPrior)
	    {
	      if(not v->priorIsFitting(prior, traln))
		{
		  tout  << "You forced prior " << &prior << " to be applied to parameter  "  << v.get() << ". This is not possible. "  << std::endl; 
		  exitFunction(-1, true);
		}

	      // tout << "setting prior " << &prior << 
		// " for param "<< v.get() << std::endl; 

	      v->setPrior( std::unique_ptr<AbstractPrior>(prior.clone()) ) ; 
	      break; 
	    }
	}

      // try to find a general prior 
      if(not foundPrior)
	{
	  for(auto iter = priors.find(cat) ; iter != priors.end() ; ++iter)
	    {
	      auto &prior = *(std::get<1>(iter->second).get()); 

	      if(not v->priorIsFitting(prior, traln))
		continue; 

	      auto &partitionsOfPrior = std::get<0>(iter->second); 
	      foundPrior = partitionsOfPrior.size() ==  0; 
	      if(foundPrior)
		{
		  // tout << "setting prior " << &prior << " for param "<< v.get() << std::endl; 
		  v->setPrior( std::unique_ptr<AbstractPrior>(prior.clone()) ) ; 
		  break; 
		}
	    }
	}
      
      if(not foundPrior)
	addStandardPrior(v.get(),traln); 
    }
}


void RunFactory::addSecondaryParameters(AbstractProposal* proposal, const std::vector<unique_ptr<AbstractParameter> > &allParameters)
{
  // get all branch length pararemeters   
  auto blParameters = std::vector<unique_ptr<AbstractParameter> >{}; 
  for(auto &v : allParameters)
    if(v->getCategory() == Category::BRANCH_LENGTHS)
      blParameters.push_back(std::unique_ptr<AbstractParameter>(v->clone())); 

  bool needsBl = false; 
  std::unordered_set<nat> myPartitions; 
  for(auto &v: proposal->getPrimaryParameterView())
    {
      needsBl |=  ( v->getCategory() == Category::FREQUENCIES 
		   || v->getCategory() == Category::SUBSTITUTION_RATES
		   || v->getCategory() == Category::TOPOLOGY
		   || v->getCategory() == Category::AA_MODEL
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



std::tuple<std::vector<std::unique_ptr<AbstractProposal>>,std::vector<ProposalSet> > 
 RunFactory::produceProposals(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, 
			      std::vector<std::unique_ptr<AbstractParameter> > &params, const TreeAln &traln, bool componentWiseMH)
{
  auto proposals = std::vector<std::unique_ptr<AbstractProposal> >{} ; 
  auto resultPropSet = std::vector<ProposalSet >{} ; 
  addPriorsToParameters(traln, priorInfo, params);

  auto reg = ProposalRegistry{}; 

  // instantiate all proposals that integrate over one parameter 
  for(auto &v : params)
    {
      if(not v->getPrior()->needsIntegration() )
	continue;

      if( componentWiseMH && v->getCategory() != Category::TOPOLOGY )
	continue; 
      
      auto tmpResult = reg.getSingleParameterProposals(v->getCategory(), propConfig,  traln); 

      // remove proposals that are not meant for DNA/AA
      if(v->getCategory() == Category::SUBSTITUTION_RATES)
	{
	  // bool isProtPartition = traln.getPartition(v->getPartitions().at(0)).dataType == AA_DATA; 
	  bool isDNAPartition = traln.getPartition(v->getPartitions().at(0)).getDataType() == PLL_DNA_DATA; 
	  auto tmp = decltype(tmpResult){}; 
	  for(auto &elem : tmpResult )
	    {
	      if( not ( elem->isForProteinOnly() && isDNAPartition ) ) // TODO more generic!  
		tmp.push_back(std::move(elem)); 
	    }
	  tmpResult.clear(); 
	  tmpResult = std::move( tmp ) ; 
	}
      
      for(auto &p : tmpResult )
	{
	  p->addPrimaryParameter(std::unique_ptr<AbstractParameter>(v->clone()));
	  addSecondaryParameters(p.get(), params); 
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
  auto mashableParameters = std::vector<AbstractParameter*>{}; 
  for(auto &v : params)
    {
      if(not v->getPrior()->needsIntegration() ) 
	continue;

      // those are prototypes! non-owning pointers 
      switch(v->getCategory())
	{
	case Category::SUBSTITUTION_RATES: 
	case Category::FREQUENCIES:
	case Category::RATE_HETEROGENEITY: 
	case Category::AA_MODEL:	
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

  // TODO that's all a bit cumbersome, has to be re-designed a bit 
  if(componentWiseMH)
    {
      auto cat2param = std::unordered_map<Category, std::vector<AbstractParameter*> >{}; 
      for(auto &p:  mashableParameters)
	cat2param[p->getCategory()].push_back(p); 

      for(auto &elem : cat2param)
	{
	  auto proposalsForSet = reg.getSingleParameterProposals(elem.first, propConfig, traln);

	  // filter out proposals that do not fit the current data
	  // type of
	  // TODO improve the setup here 
	  if(std::get<0>(elem) == Category::SUBSTITUTION_RATES)
	    {
	      // bool isProtPartition = traln.getPartition(std::get<1>(elem).at(0)->getPartitions().at(0)).dataType == AA_DATA; 
	      bool isDNAPartition = traln.getPartition(std::get<1>(elem).at(0)->getPartitions().at(0)).getDataType() == PLL_DNA_DATA; 
	      auto tmp = decltype(proposalsForSet){}; 
	      for(auto &elem : proposalsForSet )
		{
		  if( not ( elem->isForProteinOnly() &&  isDNAPartition ) ) // TODO more generic! 
		    tmp.push_back(std::move(elem)); 
		}
	      proposalsForSet.clear(); 
	      proposalsForSet = std::move( tmp ) ; 
	    }	  
	  
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
		  addSecondaryParameters(proposalClone.get(), params); 
		  lP.push_back(std::move(proposalClone));		  
		} 
	      proposalType->setInSetExecution(true);

	      resultPropSet.emplace_back(lP[0]->getRelativeWeight(), std::move(lP)); 
	    }
	}
    }

  return std::make_tuple(std::move(proposals), resultPropSet); 
}

