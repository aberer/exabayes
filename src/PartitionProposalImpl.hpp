
template<typename FUN, typename PARAM>
PartitionProposal<FUN,PARAM>::PartitionProposal(Chain *_chain, double relativeWeight, double _param, string _name)
  :  parameter(_param)
  ,chain(_chain)  
{
  this->relativeProbability = relativeWeight;
  this->name= _name;
  this->category = PARAM::cat; 
  // ptype? 
}



template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) 
{
  model = rand.drawRandInt( traln.getNumberOfPartitions());  
  values = PARAM::getParameters(traln, model); 
  vector<double> proposedValues =  FUN::getNewValues(values, parameter, rand, hastings); 
  PARAM::setParameters(traln, model, proposedValues);
  PARAM::init(traln, model);
}

template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::evaluateProposal(TreeAln &traln, PriorManager &prior) 
{
  branch root = findRoot(chain->traln->getTr());  
  nodeptr p = findNodeFromBranch(chain->traln->getTr(), root); 
  evaluateOnePartition(chain, p, TRUE, model); 
}
  

template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::resetState(TreeAln &traln, PriorManager &prior) 
{
  PARAM::setParameters(traln, model, values);
  PARAM::init(traln,model);
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::autotune() 
{
  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, FALSE);
  
#ifdef DEBUG_PRINT_TUNE_INFO
  cout << pf->name << ": with ratio " << ctr->getRatioInLastInterval() << ": "<< ((newParam < *parameter ) ? "reducing" : "increasing") <<  "\t" << *parameter << "," << newParam << endl; 
#endif
  
  parameter = newParam; 
  sctr.nextBatch();
}



