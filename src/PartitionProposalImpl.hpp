
template<typename FUN, typename PARAM>
PartitionProposal<FUN,PARAM>::PartitionProposal(double relativeWeight, double _param, string _name)
  :  parameter(_param)
{
  this->relativeProbability = relativeWeight;
  this->name= _name;
  this->category = PARAM::cat; 
  // ptype? 
}



template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  model = rand.drawRandInt( traln.getNumberOfPartitions());  
  values = PARAM::getParameters(traln, model); 
  vector<double> proposedValues =  FUN::getNewValues(values, parameter, rand, hastings); 
  PARAM::setParameters(traln, model, proposedValues);
  PARAM::init(traln, model);
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
  PARAM::setParameters(traln, model, values);
  PARAM::init(traln,model);
}



template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::autotune() 
{
  if(not FUN::tune)
    return; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, not FUN::tuneup);
  
#ifdef DEBUG_PRINT_TUNE_INFO
  cout << pf->name << ": with ratio " << ctr->getRatioInLastInterval() << ": "<< ((newParam < *parameter ) ? "reducing" : "increasing") <<  "\t" << *parameter << "," << newParam << endl; 
#endif
  
  parameter = newParam; 
  sctr.nextBatch();
}


template<typename FUN, typename PARAM> 
PartitionProposal<FUN,PARAM>* PartitionProposal<FUN,PARAM>::clone() const
{
  return new PartitionProposal<FUN,PARAM>( relativeProbability, parameter, name);
}







