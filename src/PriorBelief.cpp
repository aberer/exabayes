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

  cout << "initial prior is " << lnPrior << endl;  
}


void PriorBelief::accountForFracChange(const TreeAln &traln, int model, const vector<double> &oldFc, const vector<double> &newFcs )  
{
  assert(0); 
}


double PriorBelief::scoreEverything(const TreeAln &traln) const 
{
  double result = 0; 
 
  for(auto v : variables)
    {
      switch(v.getCategory())			// TODO => category object 
	{	  
	case TOPOLOGY: 	  
	  // result += v.getPrior()->getLogProb({});
	  result += 0; 	// well...
	  break; 
	case BRANCH_LENGTHS: 
	  {
	    assert(traln.getNumBranches() == 1); 
	    vector<branch> bs; 
	    extractBranches(traln,bs);
	    auto pr = v.getPrior();	    
	    for(auto b : bs)
	      result += pr->getLogProb( { branchLengthToReal(traln.getTr(),b.length[0]) } );
	  }
	  break; 
	case FREQUENCIES: 
	  {
	    auto freqs = traln.getFrequencies(v.getPartitions()[0]); 
	    result += v.getPrior()->getLogProb( freqs) ; 
	  } 
	  break; 
	case SUBSTITUTION_RATES: 
	  {
	    auto revMat = traln.getRevMat(v.getPartitions()[0]); 
	    result += v.getPrior()->getLogProb(revMat); 
	  }
	  break; 
	case RATE_HETEROGENEITY: 
	  {
	    double alpha = traln.getAlpha(v.getPartitions()[0]) ;
	    result += v.getPrior()->getLogProb({ alpha }); 
	  }
	  break; 
	case AA_MODEL: 
	  assert(0); 
	  break; 
	default : assert(0); 
	}
    }

  return result; 
} 


void PriorBelief::verifyPrior(const TreeAln &traln) const 
{
  double verified = scoreEverything(traln); 
  assert(verified == lnPrior); 
}


void PriorBelief::updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, shared_ptr<AbstractPrior> brPr) 
{
  if(typeid(brPr) == typeid(shared_ptr<ExponentialPrior>))
    lnPrior += (newInternalZ - oldInternalZ) *  traln.getTr()->fracchange; 
  else 
    lnPrior += brPr->getLogProb({newInternalZ}) - brPr->getLogProb({oldInternalZ}) ; 
} 


