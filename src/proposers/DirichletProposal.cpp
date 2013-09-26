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

  auto newValues = std::vector<double> (); 
  auto scalingFactor = parameter * oldValues.size(); 
  
  auto scaledOld = oldValues; 
  for(auto &v : scaledOld)
    v *= scalingFactor; 
  
  newValues = rand.drawRandDirichlet(scaledOld); 

  if( minMaxIsRelative)		// for revmat 
    {
      for(auto &v : newValues)
      	v /= *(newValues.rbegin()); 
      
      BoundsChecker::correctRevMat(newValues);       

      double sum = std::accumulate(newValues.begin(), newValues.end(), 0.); 
      for(auto &v : newValues)
      	v /= sum; 
    }
  else 				// for freuqencies (or everything else, that sums up to 1)
    {      
      correctAbsoluteRates(newValues); 
    }  

  auto scaledNew = newValues; 
  for(auto &v : scaledNew)
    v *= scalingFactor; 

  auto backP = Density::lnDirichlet(oldValues, scaledNew); 
  auto forP = Density::lnDirichlet(newValues, scaledOld); 
  
  AbstractProposal::updateHastingsLog(hastings, backP - forP   , "dirichlet"); 
    
  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 
  return newValues; 
}
