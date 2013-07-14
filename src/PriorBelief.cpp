#include "PriorBelief.hpp"
#include "Branch.hpp"
#include "GlobalVariables.hpp"
#include "Priors.hpp"

#include "Category.hpp"



PriorBelief::PriorBelief()
  :lnPrior(0)
  ,lnPriorRatio(0)
  , wasInitialized(false)
{
}


void PriorBelief::initialize(const TreeAln &traln, const std::vector<AbstractParameter*> &variables)
{
  lnPrior = scoreEverything(traln, variables); 
  lnPriorRatio = 0; 
  wasInitialized = true; 
}


void PriorBelief::accountForFracChange(const TreeAln &traln, const std::vector<double> &oldFc, const std::vector<double> &newFcs, const std::vector<AbstractPrior* > &blPriors)  
{
  assert(wasInitialized); 
  assert(blPriors.size() == 1  &&  dynamic_cast<ExponentialPrior*>(blPriors[0]) != nullptr) ; 

  // TODO this is horrible, but let's go with that for now. 

  double lambda = dynamic_cast<ExponentialPrior*> (blPriors[0])->getLamda(); 

#ifdef EFFICIENT
  // TODO investigate on more efficient tree length 
  assert(0); 
#endif

  assert(oldFc.size() == 1 && newFcs.size() == 1 );  

  double blInfluence = 1; 
  std::vector<Branch> branches = traln.extractBranches() ;  
  for(auto &b : branches)
    blInfluence *= b.getLength(); 

  lnPriorRatio += (newFcs[0] - oldFc[0]) * lambda * log(blInfluence);
}


double PriorBelief::scoreEverything(const TreeAln &traln, std::vector<AbstractParameter*> variables) const 
{
  double result = 0; 

  // for(auto &v = variables.begin(); v < variables.end(); ++v)  
       for( const auto v : variables)
    {
      double partialResult = 0; 

      switch(v->getCategory())			// TODO => category object 
	{	  
	case Category::TOPOLOGY: 	  
	  // partialResult = v.getPrior()->getLogProb({});
	  partialResult = 0; 	// well...
	  break; 
	case Category::BRANCH_LENGTHS: 
	  {
	    assert(traln.getNumBranches() == 1); 
	    std::vector<Branch> bs = traln.extractBranches(); 
	    auto pr = v->getPrior();	    
	    for(auto b : bs)
	      partialResult += pr->getLogProb( {  b.getInterpretedLength(traln)} );
	  }
	  break; 
	case Category::FREQUENCIES: 
	  {
	    auto freqs = traln.getFrequencies(v->getPartitions()[0]); 
	    partialResult = v->getPrior()->getLogProb( freqs) ; 
	  } 
	  break; 
	case Category::SUBSTITUTION_RATES: 
	  {
	    auto revMat = traln.getRevMat(v->getPartitions()[0]); 
	    partialResult = v->getPrior()->getLogProb(revMat); 
	  }
	  break; 
	case Category::RATE_HETEROGENEITY: 
	  {
	    double alpha = traln.getAlpha(v->getPartitions()[0]) ;
	    partialResult = v->getPrior()->getLogProb({ alpha }); 
	  }
	  break; 
	case Category::AA_MODEL: 
	  assert(0); 
	  break; 
	default : assert(0); 
	}
      
      // cout << "pr(" <<  v << ") = "<< partialResult << endl ; 
      result += partialResult; 
    }

  return result; 
} 


void PriorBelief::verifyPrior(const TreeAln &traln, std::vector<AbstractParameter*> variables) const 
{
  assert(lnPriorRatio == 0); 
  double verified = scoreEverything(traln, variables); 
  if ( fabs(verified -  lnPrior ) >= ACCEPTED_LNPR_EPS)
    {
      std::cerr << std::setprecision(10) << "ln prior was " << lnPrior << " while it should be " << verified << std::endl; 
      assert(fabs(verified -  lnPrior ) < ACCEPTED_LNPR_EPS); 
    }
}


void PriorBelief::updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, AbstractPrior* brPr) 
{
  assert(wasInitialized); 
  if(dynamic_cast<ExponentialPrior*> (brPr) != nullptr ) 
    {
      auto casted = dynamic_cast<ExponentialPrior*> (brPr); 
      lnPriorRatio += (log(newInternalZ /   oldInternalZ) ) * traln.getTr()->fracchange * casted->getLamda(); 
    }
  else 
    {
      assert(0);		// very artificial
      lnPriorRatio += brPr->getLogProb({newInternalZ}) - brPr->getLogProb({oldInternalZ}) ; 
    }
}



