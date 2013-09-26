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
      correctAbsoluteRates(oldValues);
    }


  double sum = std::accumulate(oldValues.begin(), oldValues.end(), 0.); 
  if( fabs(sum - 1.0 ) > 1e-6 )
    {
      std::cerr << "Danger: while proposing values,  sum was " << sum << ". values: "   << oldValues  << std::endl; 
      assert(0); 
    }

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
