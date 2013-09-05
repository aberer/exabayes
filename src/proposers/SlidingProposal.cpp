#include "SlidingProposal.hpp"
#include "BoundsChecker.hpp"

SlidingProposal::SlidingProposal(double minVal, double maxVal, bool _minMaxIsRelative)
  : minMaxIsRelative(_minMaxIsRelative)
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

  // tout << "proposing for " << oldValues << std::endl; 

  double oldVal = oldValues.at(posA) ; 
  double newVal = rand.drawFromSlidingWindow(oldVal, parameter); 

  if( newVal <= 0)
    {
      // append to the upper end of the window 
      newVal = oldVal + parameter / 2  - newVal; 
      assert(newVal >= 0 && newVal < 1); 
    }
  else if (newVal >= 1)
    {
      // append to the lower end of the window 
      newVal = oldVal - parameter / 2 - newVal + 1 ; 
      assert(newVal >= 0 && newVal < 1); 
    }


  oldValues[posA] = newVal; 


  if(minMaxIsRelative)
    {
      for(auto &v : oldValues)
	v /= *(oldValues.rbegin()); 

      BoundsChecker::correctRevMat(oldValues); 

      double sum = std::accumulate(oldValues.begin(), oldValues.end(), 0.); 
      for(auto &v : oldValues) 
	v /= sum; 
    }
  else 
    {
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
    }

  assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0.)  - 1.0 ) < 1e-6 ); 

  return oldValues; 
}



std::vector<double> SlidingProposal::proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
  // tout << "slider" << std::endl; 


  if(oldValues.size() == 1 )
    return {proposeOneValue(oldValues[0], parameter,rand, hastings)}; 
  else
    return proposeRelativeMany(oldValues, parameter, rand, hastings);
}
