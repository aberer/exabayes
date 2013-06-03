#include "PriorBelief.hpp"
#include "branch.h"
#include "GlobalVariables.hpp"
#include "Priors.hpp"


PriorBelief::PriorBelief(const TreeAln &traln, const vector<RandomVariable> &_variables)
  :lnPrior(0)
  , lnPriorRatio(0)
  , variables(_variables)
{
  lnPrior = scoreEverything(traln);
}


void PriorBelief::accountForFracChange(const TreeAln &traln, int model, const vector<double> &oldFc, const vector<double> &newFcs, double lambda )  
{
  assert(oldFc.size() == 1 && newFcs.size() == 1 );  
  lnPriorRatio += (newFcs[0] - oldFc[0]) * lambda * traln.getTreeLength();
}


double PriorBelief::scoreEverything(const TreeAln &traln) const 
{
  double result = 0; 

  for(auto v : variables)
    {
      double partialResult = 0; 

      switch(v.getCategory())			// TODO => category object 
	{	  
	case TOPOLOGY: 	  
	  // partialResult = v.getPrior()->getLogProb({});
	  partialResult = 0; 	// well...
	  break; 
	case BRANCH_LENGTHS: 
	  {
	    assert(traln.getNumBranches() == 1); 
	    vector<branch> bs; 
	    extractBranches(traln,bs);
	    auto pr = v.getPrior();	    
	    for(auto b : bs)
	      partialResult += pr->getLogProb( { branchLengthToReal(traln.getTr(),b.length[0]) } );
	  }
	  break; 
	case FREQUENCIES: 
	  {
	    auto freqs = traln.getFrequencies(v.getPartitions()[0]); 
	    partialResult = v.getPrior()->getLogProb( freqs) ; 
	  } 
	  break; 
	case SUBSTITUTION_RATES: 
	  {
	    auto revMat = traln.getRevMat(v.getPartitions()[0]); 
	    partialResult = v.getPrior()->getLogProb(revMat); 
	  }
	  break; 
	case RATE_HETEROGENEITY: 
	  {
	    double alpha = traln.getAlpha(v.getPartitions()[0]) ;
	    partialResult = v.getPrior()->getLogProb({ alpha }); 
	  }
	  break; 
	case AA_MODEL: 
	  assert(0); 
	  break; 
	default : assert(0); 
	}
      
      // cout << "pr(" <<  v << ") = "<< partialResult << endl ; 
      result += partialResult; 
    }

  return result; 
} 


void PriorBelief::verifyPrior(const TreeAln &traln) const 
{
  assert(lnPriorRatio == 0); 
  double verified = scoreEverything(traln); 
  if ( fabs(verified -  lnPrior ) >= 1e-6)
    {
      cerr << "ln prior was " << lnPrior << " while it should be " << verified << endl; 
      assert(fabs(verified -  lnPrior ) < 1e-6); 
    }
}

shared_ptr<AbstractPrior> PriorBelief::getBranchLengthPrior() const
{
  for(auto v : variables)
    {
      if(v.getCategory() == BRANCH_LENGTHS)
	{
	  // cout << "prior is " << v.getPrior() << endl; 
	  return v.getPrior(); 
	}
    }

  assert(0); 
  return shared_ptr<AbstractPrior>(); 
}



void PriorBelief::updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, shared_ptr<AbstractPrior> brPr) 
{
  if(dynamic_cast<ExponentialPrior*> (brPr.get()) != nullptr ) 
    {
      auto casted = dynamic_cast<ExponentialPrior*> (brPr.get()); 
      lnPriorRatio += (log(newInternalZ /   oldInternalZ) ) * traln.getTr()->fracchange * casted->getLamda(); 
    }
  else 
    {
      assert(0);		// very artificial
      lnPriorRatio += brPr->getLogProb({newInternalZ}) - brPr->getLogProb({oldInternalZ}) ; 
    }
} 


