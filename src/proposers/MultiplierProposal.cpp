#include "MultiplierProposal.hpp"


MultiplierProposal::MultiplierProposal(double minVal, double maxVal)
{
  this->minVal = minVal; 
  this->maxVal = maxVal; 
  tune = true; 
  tuneup = true; 
}


std::vector<double> MultiplierProposal::proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{    
  double newVal = 0, position = 0, multiplier = 0; 

  position = rand.drawIntegerOpen(oldValues.size()); 
  multiplier =  rand.drawMultiplier( parameter);     
  newVal = oldValues[position] * multiplier; 

  // TODO allowed? 
  if(newVal < minVal)
    newVal = minVal; 
  else if(maxVal < newVal)
    newVal = maxVal; 

  AbstractProposal::updateHastingsLog(hastings, log(multiplier), "multiplier"); 
  oldValues[position] = newVal; 
  return oldValues;
}
