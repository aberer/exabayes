#include "DirichletProposal.hpp"
#include "BoundsChecker.hpp" 

DirichletProposal::DirichletProposal( double minVal, double maxVal, bool _minMaxIsRelative ) 
  : minMaxIsRelative(_minMaxIsRelative)
{
  this->minVal = minVal; 
  this->maxVal = maxVal; 
  tune = true ; 
  tuneup = false; 
} 


std::vector<double> DirichletProposal::proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 

  // tout << "dirichlet" << std::endl; 

  auto newValues = std::vector<double> (); 
    
  auto scaledOld = oldValues; 
  for_each(scaledOld.begin(), scaledOld.end(), [&](double &d) { return d *= parameter * oldValues.size() ; }) ; 
  newValues = rand.drawRandDirichlet(scaledOld); 

  if( minMaxIsRelative)		// for revmat 
    {
      // tout << MAX_SCI_PRECISION << "initial "  << newValues << std::endl; 
      for(auto &v : newValues)
      	v /= *(newValues.rbegin()); 
      
      BoundsChecker::correctRevMat(newValues);       

      // tout << MAX_SCI_PRECISION << "after correction "  << newValues << std::endl; 
      
      double sum = std::accumulate(newValues.begin(), newValues.end(), 0.); 
      for(auto &v : newValues)
      	v /= sum; 
      // tout << MAX_SCI_PRECISION << "after check "  << newValues << std::endl; 
    }
  else 				// for freuqencies (or everything else, that sums up to 1)
    {      
      // correct for problematic values 
      int numHigh = 0,
	numLow = 0; 
      double normer = 0; 
      for(auto &v : newValues) 
	{
	  if(v < minVal)
	    {
	      v = minVal; 
	      ++numLow ; 
	    }
	  else if(maxVal < v)
	    {
	      v = maxVal; 
	      ++numHigh; 
	    }
	  else 
	    normer += v ; 
	}

      normer =   (1. -( numHigh *maxVal + numLow * minVal))  / normer; 
    
      for(auto &v : newValues)
	{
	  if( v != minVal && v != maxVal)
	    v *= normer; 
	}
    }  

  AbstractProposal::updateHastings(hastings, densityDirichlet(oldValues, newValues) / densityDirichlet(newValues,oldValues) , "dirichlet"); 
    
  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 
  return newValues; 
}
