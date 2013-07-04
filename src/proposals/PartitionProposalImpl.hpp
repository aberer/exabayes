
template<typename FUN, typename PARAM >
PartitionProposal<FUN,PARAM>::PartitionProposal( double _param, string _name)
  :  parameter(_param)
{  
  this->name= _name;
  this->category = PARAM::cat; 
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  auto oneModel = primVar[0]->getPartitions()[0];
  assert(primVar.size() == 1 ) ; // <= TODO 

  values = PARAM::getParameters(traln, oneModel); 
  double oldFracChange = traln.getTr()->fracchange; 
  vector<double> proposedValues; 

  proposedValues = FUN::getNewValues(values, parameter, rand, hastings); 

  assert(proposedValues.size() == values.size()); 
  assert(traln.getNumBranches() == 1 ); 

  for(auto v : primVar[0]->getPartitions())
    PARAM::setParameters(traln, v, proposedValues);  

  double newFracChange = traln.getTr()->fracchange;   
  assert(primVar.size() != 0);   

  if(PARAM::modifiesBL)
    {
      vector<AbstractPrior* > blPriors; 
      for(auto v : secVar)
	blPriors.push_back(v->getPrior()) ; 
      prior.accountForFracChange(traln, {oldFracChange}, {newFracChange} , blPriors);
    }
 
  auto thePrior = primVar[0]->getPrior(); 
  prior.addToRatio( thePrior->getLogProb(proposedValues) -  thePrior->getLogProb(values) ) ;
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior) 
{
  assert(primVar.size() == 1 ); 
  evaluator.evaluatePartitions(traln, primVar[0]->getPartitions() ); 
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::resetState(TreeAln &traln, PriorBelief &prior) 
{
  assert(primVar.size() == 1 ); 

  for(auto elem : primVar[0]->getPartitions())
    {
      PARAM::setParameters(traln, elem, values);
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
