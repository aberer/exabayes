#include "PriorBelief.hpp"

#include <cmath>

#include "Branch.hpp"
#include "GlobalVariables.hpp"
#include "priors/AbstractPrior.hpp"
#include "priors/UniformPrior.hpp"
#include "priors/ExponentialPrior.hpp"
#include "Category.hpp"

#include "TreePrinter.hpp"


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


void PriorBelief::accountForFracChange(TreeAln &traln, const std::vector<double> &oldFcs, const std::vector<double> &newFcs, 
				       const std::vector<AbstractParameter*> &affectedBlParams )  
{
  assert(wasInitialized); 

  nat ctr = 0; 
  for(auto &param : affectedBlParams)
    {
      lnPriorRatio += param->getPrior()->accountForMeanSubstChange(traln,  param, oldFcs.at(ctr), newFcs.at(ctr));
      ++ctr; 
    }
}


double PriorBelief::scoreEverything(const TreeAln &traln, const std::vector<AbstractParameter*> &parameters) const 
{
  double result = 0; 

  for( auto& v : parameters)
    {
      double partialResult = 0; 

      switch(v->getCategory()) 
	{	  
	case Category::TOPOLOGY: 	  
	  partialResult = 0; 	// well...
	  break; 
	case Category::BRANCH_LENGTHS: 
	  {
	    auto bs = traln.extractBranches(v); 
	    auto pr = v->getPrior();	    
	    for(auto &b : bs)	      
	      partialResult += pr->getLogProb( ParameterContent {  { b.getInterpretedLength(traln,v) }  }); 
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
	    auto revMat = traln.getRevMat(v->getPartitions()[0] , false); 
	    partialResult = v->getPrior()->getLogProb(revMat); 
	  }
	  break; 
	case Category::RATE_HETEROGENEITY: 
	  {
	    double alpha = traln.getAlpha(v->getPartitions()[0]) ;
	    partialResult = v->getPrior()->getLogProb( ParameterContent{{ alpha} }); 
	  }
	  break; 
	case Category::AA_MODEL: 
	  {
	    auto p = v->getPartitions()[0]; 
	    auto model = traln.getProteinModel(p); 
	    auto content =  ParameterContent(); 
	    content.protModel.push_back(model); 
	    partialResult = v->getPrior()->getLogProb(content);
	  }
	  break; 
	default : assert(0); 
	}
      
      // std::cout << "pr(" <<  v << ") = "<< partialResult << std::endl ; 
      result += partialResult; 
    }

  
  // the lnPriorRatio may be infinite. But never the assumed value 
  assert(not std::isinf(result) && not std::isnan(result)); 
  
  return result; 
} 


void PriorBelief::verifyPrior(const TreeAln &traln, std::vector<AbstractParameter*> parameters) const 
{
  assert(lnPriorRatio == 0); 
  double verified = scoreEverything(traln, parameters); 
  if ( fabs(verified - lnPrior) >= ACCEPTED_LNPR_EPS)
    {
      tout << MAX_SCI_PRECISION << "ln prior was " << lnPrior << " while it should be " << verified << std::endl; 
      tout << "difference: "<< fabs(verified - lnPrior) << std::endl; 
      assert(fabs(verified -  lnPrior ) < ACCEPTED_LNPR_EPS); 
    }
}

