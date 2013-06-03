
template<typename FUN, typename PARAM>
PartitionProposal<FUN,PARAM>::PartitionProposal( double _param, string _name)
  :  parameter(_param)
{  
  this->name= _name;
  this->category = PARAM::cat; 
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  auto oneModel  = randomVariables[0].getPartitions()[0];

  cout <<  "we have " << randomVariables.size() << " random variables " << endl; 
  for(auto v : randomVariables)
    cout <<  v << endl; 

  assert(randomVariables.size() ==1 ) ; // <= TODO 

  values = PARAM::getParameters(traln, oneModel); 
  vector<double> proposedValues =  FUN::getNewValues(values, parameter, rand, hastings); 
  assert(proposedValues.size() == values.size()); 

  assert(traln.getNumBranches() == 1 ); 

  double oldFracChange = traln.getTr()->fracchange; 

  for(auto v : randomVariables[0].getPartitions())
    {
      PARAM::setParameters(traln, v, proposedValues);  
      PARAM::init(traln, v);  
    }

  double newFracChange = traln.getTr()->fracchange;   

  assert(randomVariables.size() != 0);   

  auto thePrior = randomVariables[0].getPrior(); 
  prior.addToRatio( thePrior->getLogProb(proposedValues) -  thePrior->getLogProb(values) ) ;
  
  auto brPr = prior.getBranchLengthPrior();
  if(PARAM::needsFcUpdate )
    {
      if( dynamic_cast<ExponentialPrior*>(brPr.get()) != nullptr)
	{	  
	  double lambda = dynamic_cast<ExponentialPrior*>(brPr.get())->getLamda();
	  prior.accountForFracChange(traln, oneModel, {oldFracChange}, {newFracChange}, lambda);
	}
      else if (dynamic_cast<UniformPrior*>(brPr.get()) != nullptr) 
	{
	  // uniform prior: all branches need to be checked =/ 
	  assert(0); 
	}
      else 
	assert(0); 
    }
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  // branch root = findRoot(traln.getTr());  

  // nodeptr p = findNodeFromBranch(traln.getTr(), root); 

  assert(randomVariables.size() == 1 ); 

#ifdef EFFICIENT
  // dammit, we need to do it with the root...
  assert(0); 
#endif
  // evaluatePartitions(traln, traln.getTr()->start,   randomVariables[0].getPartitions()) ; 
  evaluateGenericWrapper(traln, traln.getTr()->start, TRUE); 

}

template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::resetState(TreeAln &traln, PriorBelief &prior) 
{
  assert(randomVariables.size() == 1 ); 

  for(auto elem : randomVariables[0].getPartitions())
    {
      PARAM::setParameters(traln, elem, values);
      PARAM::init(traln,elem);  
    }
}



template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::autotune() 
{
  if(not FUN::tune)
    return; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, not FUN::tuneup);
  
#ifdef DEBUG_PRINT_TUNE_INFO
  cout << name << ": with ratio " << sctr.getRatioInLastInterval() << ": "<< ((newParam < parameter ) ? "reducing" : "increasing") <<  "\t" << parameter << "," << newParam << endl; 
#endif
  
  parameter = newParam; 
  
  sctr.nextBatch();
}


template<typename FUN, typename PARAM> 
PartitionProposal<FUN,PARAM>* PartitionProposal<FUN,PARAM>::clone() const
{
  return new PartitionProposal<FUN,PARAM>(  *this );
}
