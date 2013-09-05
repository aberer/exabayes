#include "PriorBelief.hpp"

#include "Branch.hpp"
#include "GlobalVariables.hpp"
#include "priors/AbstractPrior.hpp"
#include "priors/ExponentialPrior.hpp"
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


void PriorBelief::accountForFracChange(const TreeAln &traln, const std::vector<double> &oldFc, const std::vector<double> &newFcs, 
				       const std::vector<AbstractParameter*> &affectedBlParams )  
{
  // TODO this is horrible, but let's go with that for now. 

  assert(wasInitialized); 

  for(auto &p :  affectedBlParams)
    {
      // TODO if fixed, we really have to go through all branch lengths
      // and change them      
      // => disabled currently  
      auto tmp = p->getPrior();
      assert(dynamic_cast<ExponentialPrior*>(tmp) != nullptr) ; 
    }

#ifdef EFFICIENT
  // TODO investigate on more efficient tree length 
  assert(0); 
#endif


  for(auto &param : affectedBlParams)
    {
      auto branches = traln.extractBranches(param); 

      double lambda = dynamic_cast<ExponentialPrior*> (param->getPrior())->getLamda();
      double myOld = oldFc.at(param->getIdOfMyKind()); 
      double myNew = newFcs.at(param->getIdOfMyKind()); 

      double blInfluence = 1; 
      for(auto &b : branches)
	blInfluence *= b.getLength();
	lnPriorRatio += (myNew - myOld) * lambda * log(blInfluence);
      
    }

  // auto branches = traln.extractBranches(affectedBlParams); 

  // for(auto &param : affectedBlParams)
  //   {
  //     double lambda = dynamic_cast<ExponentialPrior*> (param->getPrior())->getLamda();
  //     double myOld = oldFc.at(param->getIdOfMyKind()); 
  //     double myNew = newFcs.at(param->getIdOfMyKind()); 

  //     double blInfluence = 1; 
  //     for(auto &b : branches)
  // 	blInfluence *= b.getLength(param);

  // 	lnPriorRatio += (myNew - myOld) * lambda * log(blInfluence);
  //   }
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
	    std::vector<BranchLength> bs = traln.extractBranches(v); 
	    auto pr = v->getPrior();	    
	    for(auto &b : bs)	      
	      {
		partialResult += pr->getLogProb( {  b.getInterpretedLength(traln,v) }); 
	      }
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
      
      // std::cout << "pr(" <<  v << ") = "<< partialResult << std::endl ; 
      result += partialResult; 
    }
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


void PriorBelief::updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,
					  double newInternalZ, const AbstractParameter* param ) 
{
  auto prior = param->getPrior();
  assert(wasInitialized); 
  if(dynamic_cast<ExponentialPrior*> (prior) != nullptr ) 
    {
      auto casted = dynamic_cast<ExponentialPrior*> (prior); 
      lnPriorRatio += 
	(log(newInternalZ /   oldInternalZ) ) 
	* traln.getMeanSubstitutionRate(param->getPartitions()) * casted->getLamda(); 
    }
  else 
    {
      assert(0);		// very artificial
      lnPriorRatio += prior->getLogProb({newInternalZ}) - prior->getLogProb({oldInternalZ}) ; 
    }
}



