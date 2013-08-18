#include "AlignmentProposal.hpp"

#define HEAT 0.5

#define PRINT_INFO 

#define WEIGHT_EPSILON  0.000001 


AlignmentProposal::AlignmentProposal(Category cat, std::string name, double parameter, 
				     nat numParamExpected, AbstractProposer* proposalPrototype ) 
  : AbstractProposal()
  , partitionParameter(numParamExpected, parameter)    
{
  this->name = name; 
  this->category = cat; 
  for(nat i = 0; i < numParamExpected; ++i)
    partitionProposer.emplace_back(proposalPrototype->clone());  
}

AlignmentProposal::AlignmentProposal(const AlignmentProposal& rhs )  
  : AbstractProposal(rhs)
  , partitionParameter(rhs.partitionParameter)
{
  for(auto &p : rhs.partitionProposer)
    partitionProposer.emplace_back(p->clone());   
}


AlignmentProposal& AlignmentProposal::operator=(AlignmentProposal rhs)
{
  std::swap(*this, rhs); 
  return *this; 
}


void AlignmentProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{  
  assert(partitionProposer.size()== primaryParameters.size() 
	 && primaryParameters.size() == partitionParameter.size() ); 

  double oldFracChange = traln.getTr()->fracchange; 

  savedParams.resize(primaryParameters.size()); 
  std::vector<double> partitionLnls = traln.getPartitionLnls();
  // TODO prior ratios   
  
  std::vector<nat> affectedPartitions; 
  std::vector<double> priorRatioPerParameter; 
  std::vector<double> oldParameterLnls; 

  // propose new values for each of our partititons 
  for(nat i = 0; i < primaryParameters.size() ;++i)
    {
      auto &parameter = primaryParameters.at(i); 
      auto oldValues = parameter->extractParameter(traln); 
      savedParams[i] = oldValues; 
      auto &proposer = partitionProposer.at(i); 

      auto proposedValues = proposer->proposeValues(oldValues.values, partitionParameter.at(i), rand, hastings); 
      ParameterContent  newParameter; 
      newParameter.values = proposedValues; 
      parameter->applyParameter(traln, newParameter);       
      
      // treat the prior  
      double localPrRatio = parameter->getPrior( )->getLogProb(proposedValues) - parameter->getPrior()->getLogProb(oldValues.values); 
      priorRatioPerParameter.push_back(localPrRatio); 

      auto myPartitions = parameter->getPartitions(); 
      affectedPartitions.insert(affectedPartitions.end(), myPartitions.begin(), myPartitions.end()); 
      
      // sum up old lnls  
      double tmp = 0;       
      for(auto &p : parameter->getPartitions())
	tmp += partitionLnls[p]; 
      oldParameterLnls.push_back(tmp); 
    }

  // just a quick check, nothing should have changed a partition twice  
  std::unordered_set<nat> partitionsSeen; 
  for(auto &p : affectedPartitions)
    {
      assert(partitionsSeen.find(p) == partitionsSeen.end()) ; 
      partitionsSeen.insert(p); 
    }

  // the expensive eval-call   
  eval.evaluatePartitions(traln, affectedPartitions, true);

  // get the new lnls per parameter 
  auto newLnls = traln.getPartitionLnls(); 
  std::vector<double> newParameterLnls; 
  for(nat i = 0; i < primaryParameters.size() ;++i)
    {
      auto &parameter = primaryParameters[i]; 
      double tmp = 0; 
      for(auto &p : parameter->getPartitions())
	tmp += newLnls[p]; 
      newParameterLnls.push_back(tmp); 
    }

  // propose new parameters to the partitions according to the
  // posterior probability of the new parameters
  std::vector<nat> partitionsToReset; 
  for(nat i = 0; i < primaryParameters.size(); ++i )
    {
      auto &parameter  = primaryParameters[i]; 
      assert(priorRatioPerParameter[i] == 0.); 

      double deltaP =  exp(( oldParameterLnls[i] - newParameterLnls[i] ) - priorRatioPerParameter[i]); 

      double accProp = 1. / (1. + deltaP); // the default acceptance probability 
      
      accProp = ( accProp + WEIGHT_EPSILON ) / (1. + 2 * WEIGHT_EPSILON ) ; 

      accProp = pow(accProp,HEAT); // heat the acc prop a bit 
      
#ifdef PRINT_INFO
      tout << "for partition "<< i 
	   << "\taccProp=" << accProp
	   << "\toldLnl="             << oldParameterLnls[i]
	   << "\tnewParameterLnls=" << newParameterLnls[i]
	   << "\tpriorRatio=" << priorRatioPerParameter[i] << std::endl; 
#endif

      
      if(rand.drawRandDouble01() < accProp)
	{
#ifdef PRINT_INFO
	  tout << "partition " << i << " ACC" << std::endl; 
#endif
	  
	  double hastHere = exp(log(1-accProp)  - log( accProp)); 

	  updateHastings(hastings, hastHere, "AlnPropo" );
	  prior.addToRatio(priorRatioPerParameter[i]);
	}
      else 
	{
#ifdef PRINT_INFO
	  tout << "partition " << i << " REJ" << std::endl; 
#endif

	  double hastHere = exp(log(accProp) - log(1-accProp)) ;
	  
	  updateHastings(hastings, hastHere , "AlnPropo" ); 
	  parameter->applyParameter(traln, savedParams[i]); 
	  auto myPartitions = parameter->getPartitions(); 
	  partitionsToReset.insert(partitionsToReset.end(), myPartitions.begin(), myPartitions.end()); 	    
	}      
    }

  // reset all partitions for which we do not propose a new setting
  eval.resetSomePartitionsToImprinted(traln, partitionsToReset);
  std::vector<double> pLnlsNow = traln.getPartitionLnls();   
  for(auto &p : partitionsToReset)
    pLnlsNow.at(p) = partitionLnls.at(p);
  traln.setPartitionLnls(pLnlsNow); 


  // handle fracchange matters 
  // will be more complicated once we have per partition branch lengths 
  assert(traln.getNumBranches() == 1 ); 

  double newFracChange = traln.getTr()->fracchange; 
  if(oldFracChange != newFracChange)
    {
      std::vector<AbstractPrior*> blPriors; 
      for(auto &v : secondaryParameters)
	{
	  assert(v->getCategory() == Category::BRANCH_LENGTHS); 
	  blPriors.push_back(v->getPrior()); 
	}

      // meh 
#if 1 
      assert(0);
#else 
      // prior.accountForFracChange(traln, {oldFracChange}, {newFracChange}, blPriors); 
#endif
      updateHastings(hastings, pow(newFracChange / oldFracChange, traln.getNumberOfBranches()), name); 
    }

#ifdef PRINT_INFO
  tout << "overall hastings was " << hastings << std::endl; 
  tout << "================" << std::endl; 
#endif
  // assert(0); 
}



void AlignmentProposal::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) 
{
  // this is a NOOP  

  // applyToState must ensure that everything is okay

  auto partitionLnls = traln.getPartitionLnls();
  double overall = std::accumulate(partitionLnls.begin(), partitionLnls.end(), 0.); 
  traln.getTr()->likelihood = overall; 
}


void AlignmentProposal::resetState(TreeAln &traln)
{
  nat ctr  = 0; 
  for(auto &p : primaryParameters)
    {
      p->applyParameter(traln, savedParams[ctr]);
      ++ctr; 
    }
} 


void AlignmentProposal::autotune() 
{
  // assert(0); 
  // tout << "TODO implement autotune "  << std::endl; 
}


AbstractProposal* AlignmentProposal::clone() const 
{
  return new AlignmentProposal(*this); 
}  


void AlignmentProposal::readFromCheckpointCore(std::istream &in)
{
  for(auto &p : partitionParameter)
    p = cRead<double>(in); 
}


void AlignmentProposal::writeToCheckpointCore(std::ostream &out) const
{
  for(auto &p : partitionParameter)
    cWrite<double>(out, p); 
} 

