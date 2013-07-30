#include "SlidingProposal.hpp"

SlidingProposal::SlidingProposal(double minVal, double maxVal)
{
  this->minVal = minVal; 
  this->maxVal = maxVal; 
  tune = true; 
  tuneup = true; 
}


double SlidingProposal::proposeOneValue(double oldVal, double parameter, Randomness &rand, double &hastings)
{
  double newVal = rand.drawFromSlidingWindow(oldVal, parameter);
  if(newVal < 0 )
    newVal = - newVal ; 
  assert(0); 			// TODO  
  return newVal; 
} 


std::vector<double> SlidingProposal::proposeRelativeMany(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
  int posA = rand.drawIntegerOpen(oldValues.size()); 

  double oldVal = oldValues.at(posA) ; 
  double newVal = rand.drawFromSlidingWindow(oldVal, parameter); 

  // bouncing back and forth: even if we are above the upper limit,
  // mirroring back could bring us below the lower limit again. So
  // leather, rinse, repeat until we have a reasonable value
  while(newVal < 0 || 1 < newVal )
    {
      if(newVal < 0)
	newVal = - newVal; 
      else if(1 < newVal)
	newVal = 1 - newVal; 
    }

  oldValues[posA] = newVal; 

  double normer = 0; 
  nat numHigh = 0,
    numLow = 0; 
  for(auto &v : oldValues)
    {
      assert(v > 0); 	    
      if( v < minVal)
	{
	  v = minVal; 
	  ++numLow; 
	}
      else if(maxVal < v)
	{
	  v = maxVal; 
	  ++numHigh;
	}
      else 
	normer += v; 
    }

  normer = (1 - (minVal * numLow  +  maxVal * numHigh)) / normer ;

  for(auto &v : oldValues)
    if(v != minVal && v != maxVal)
      v *= normer; 

  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0.)  - 1.0 ) < 1e-6 ); 

  return oldValues; 
}



std::vector<double> SlidingProposal::proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
  if(oldValues.size() == 1 )
    return {proposeOneValue(oldValues[0], parameter,rand, hastings)}; 
  else
    return proposeRelativeMany(oldValues, parameter, rand, hastings);
}
