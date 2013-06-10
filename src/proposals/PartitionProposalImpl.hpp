
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
  auto oneModel = primVar[0].getPartitions()[0];
  assert(primVar.size() == 1 ) ; // <= TODO 

  values = PARAM::getParameters(traln, oneModel); 
  double oldFracChange = traln.getTr()->fracchange; 
      
  bool blCorrectionFailed = true  ; 
  vector<double> proposedValues; 
  while( blCorrectionFailed)  
    {      
      proposedValues = FUN::getNewValues(values, parameter, rand, hastings); 

      assert(proposedValues.size() == values.size()); 
      assert(traln.getNumBranches() == 1 ); 

      for(auto v : primVar[0].getPartitions())
	{
	  PARAM::setParameters(traln, v, proposedValues);  
	  PARAM::init(traln, v);  
	}

      double newFracChange = traln.getTr()->fracchange;   
      assert(primVar.size() != 0);         

      blCorrectionFailed = false; 

      // attempt correction 
      if(PARAM::needsBLupdate)
	{
	  double updateValue = oldFracChange / newFracChange; 
	  // cout << "need to update branches with " << updateValue << endl; 
	  vector<branch> oldBranches; 
	  extractBranches(traln, oldBranches); 
	  vector<branch> newBranches(oldBranches); 

	  for(auto& b : newBranches)	    
	    {
	      double &z = b.length[0]; 	      
	      z = pow(z, updateValue); 
	      
	      blCorrectionFailed = blCorrectionFailed ||  ( z < TreeAln::zMin || TreeAln::zMax < z ) ; 
	    }

	  if(not blCorrectionFailed)
	    {
	      for(auto &b : newBranches)
		{
		  nodeptr p = findNodeFromBranch(traln.getTr(), b); 
		  traln.clipNode(p,p->back, b.length[0]); 
		}

	      storedBranches = oldBranches; 
	    } 

	}
    }

  // cout << "fracchange after app " << traln.getTr()->fracchange << endl; 
 
  auto thePrior = primVar[0].getPrior(); 
  prior.addToRatio( thePrior->getLogProb(proposedValues) -  thePrior->getLogProb(values) ) ;
}


template<typename FUN, typename PARAM>
void PartitionProposal<FUN,PARAM>::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  // branch root = findRoot(traln.getTr());  

  // nodeptr p = findNodeFromBranch(traln.getTr(), root); 

  assert(primVar.size() == 1 ); 

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
  assert(primVar.size() == 1 ); 

  for(auto elem : primVar[0].getPartitions())
    {
      PARAM::setParameters(traln, elem, values);
      PARAM::init(traln,elem);  
    }

  if(PARAM::needsBLupdate)
    {
      for( auto &b : storedBranches)
	{
	  // use memset instead? 
	  nodeptr p  = findNodeFromBranch(traln.getTr(), b); 
	  traln.clipNode(p,p->back, b.length[0]); 
	}
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
