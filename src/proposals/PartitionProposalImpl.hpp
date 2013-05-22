
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
  model = rand.drawRandInt( traln.getNumberOfPartitions());  
  values = PARAM::getParameters(traln, model); 
  vector<double> proposedValues =  FUN::getNewValues(values, parameter, rand, hastings); 
  assert(proposedValues.size() == values.size()); 

  assert(traln.getNumBranches() == 1 ); 

  double oldFracChange = traln.getTr()->fracchange; 

  PARAM::setParameters(traln, model, proposedValues);  
  PARAM::init(traln, model);  

  double newFracChange = traln.getTr()->fracchange;   

  assert(randomVariables.size() == 0);   
  auto thePrior = randomVariables[0].getPrior(); 
  prior.addToRatio(   thePrior->getLogProb(proposedValues) - thePrior->getLogProb(values)  ) ;

  if(oldFracChange != newFracChange && PARAM::needsFcUpdate)
    prior.accountForFracChange(traln, model, {oldFracChange}, {newFracChange});
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  branch root = findRoot(traln.getTr());  
  nodeptr p = findNodeFromBranch(traln.getTr(), root); 
  evaluateOnePartition(traln, p, TRUE, model); 
}

template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::resetState(TreeAln &traln, PriorBelief &prior) 
{
  vector<double> curVals =  PARAM::getParameters(traln,model);
  PARAM::setParameters(traln, model, values);
  assert(curVals.size() == values.size()); 
  PARAM::init(traln,model);  
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
  // return new PartitionProposal<FUN,PARAM>(  parameter, name);
  return new PartitionProposal<FUN,PARAM>(  *this );
}
