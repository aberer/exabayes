#include "PriorBelief.hpp"

#include <cmath>

#include "model/Branch.hpp"
#include "system/GlobalVariables.hpp"
#include "priors/AbstractPrior.hpp"
#include "priors/UniformPrior.hpp"
#include "priors/ExponentialPrior.hpp"
#include "model/Category.hpp"

#include "TreePrinter.hpp"


PriorBelief::PriorBelief()
  : _lnPrior(log_double::fromAbs(1.))
  , _lnPriorRatio(log_double::fromAbs(1.))
  , wasInitialized(false)
{
}


void PriorBelief::initialize(const TreeAln &traln, const std::vector<AbstractParameter*> &variables)
{
  _lnPrior = scoreEverything(traln, variables); 
  _lnPriorRatio = log_double::fromAbs(1.); 
  wasInitialized = true; 
}


void PriorBelief::accountForFracChange(TreeAln &traln, const std::vector<double> &oldFcs, const std::vector<double> &newFcs, 
				       const std::vector<AbstractParameter*> &affectedBlParams )  
{
  assert(wasInitialized); 

  nat ctr = 0; 
  for(auto &param : affectedBlParams)
    {
      _lnPriorRatio *= log_double::fromLog(param->getPrior()->accountForMeanSubstChange(traln,  param, oldFcs.at(ctr), newFcs.at(ctr)));
      ++ctr; 
    }
}


log_double PriorBelief::scoreEverything(const TreeAln &traln, const std::vector<AbstractParameter*> &parameters) const 
{
  log_double result = log_double::fromAbs(1.); 

  for( auto& v : parameters)
    {
      log_double partialResult = log_double::fromAbs(1.); 
      
      switch(v->getCategory()) 
	{	  
	case Category::TOPOLOGY: 	  
	  partialResult *= log_double::fromAbs(1.); 	// well...
	  break; 
	case Category::BRANCH_LENGTHS: 
	  {
	    auto bs = traln.extractBranches(v); 
	    auto pr = v->getPrior();	    
	    for(auto &b : bs)	      
	      partialResult *= pr->getLogProb( ParameterContent {  { b.getInterpretedLength(traln,v) }  }); 
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
      result *= partialResult; 
    }

  // the lnPriorRatio may be infinite. But never the assumed value   
  assert(not result.isInfinity() && not result.isNaN()); 
  
  return result; 
} 


void PriorBelief::verifyPrior(const TreeAln &traln, std::vector<AbstractParameter*> parameters) const 
{
  assert(  _lnPriorRatio.toAbs() == 1.   ); 
  
  auto verified = scoreEverything(traln, parameters); 

  if (   log_double( verified / _lnPrior).toAbs() - 1.  >=  ACCEPTED_LNPR_EPS)
    {
      tout << MAX_SCI_PRECISION << "ln prior was " << _lnPrior << " while it should be " << verified << std::endl; 
      tout << "difference: "<< fabs( verified.toAbs() - _lnPrior.toAbs()) << std::endl; 
      assert(0); 
    }
}

