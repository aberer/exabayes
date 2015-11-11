
#include "DirichletProposal.hpp"


DirichletProposal::DirichletProposal( double minVal, double maxVal) 
{
  this->minVal = minVal; 
  this->maxVal = maxVal; 
  tune = true ; 
  tuneup = false; 
} 



std::vector<double> DirichletProposal::proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 

  std::vector<double> newValues; 
    
  std::vector<double> scaledOld = oldValues; 
  for_each(scaledOld.begin(), scaledOld.end(), [&](double &d) { return d *= parameter * oldValues.size() ; }) ; 
  newValues = rand.drawRandDirichlet(scaledOld); 
    
  // newValues = rand.drawDirichletExpected( oldValues, parameter * oldValues.size() );

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

  AbstractProposal::updateHastings(hastings, densityDirichlet(oldValues, newValues) / densityDirichlet(newValues,oldValues) , "dirichlet"); 
    
  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 
  return newValues; 
}
