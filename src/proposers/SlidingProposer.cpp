#include "SlidingProposer.hpp"
#include "RateHelper.hpp"
#include "BoundsChecker.hpp"
#include <numeric>

SlidingProposer::SlidingProposer(double minVal, double maxVal, bool _minMaxIsRelative)
  : AbstractProposer{true, true, minVal, maxVal}
  , minMaxIsRelative{_minMaxIsRelative}
{
}


double SlidingProposer::proposeOneValue(double oldVal, double parameter, Randomness &rand, double &hastings)
{
  double newVal = rand.drawFromSlidingWindow(oldVal, parameter);
  if(newVal < 0 )
    newVal = - newVal ; 

  // todo currently only asserting, that it is not used 
  assert(0);
  return newVal; 
} 


std::vector<double> SlidingProposer::proposeRelativeMany(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
// #if 0 
//   int posA = rand.drawIntegerOpen(oldValues.size()); 

//   double oldVal = oldValues.at(posA) ; 
//   double newVal = rand.drawFromSlidingWindow(oldVal, parameter); 

//   if( newVal <= 0)
//     {
//       // append to the upper end of the window 
//       newVal = oldVal + parameter / 2  - newVal; 
//       assert(newVal >= 0 && newVal < 1); 
//     }
//   else if (newVal >= 1)
//     {
//       // append to the lower end of the window 
//       newVal = oldVal - parameter / 2 - newVal + 1 ; 
//       assert(newVal >= 0 && newVal < 1); 
//     }

//   oldValues[posA] = newVal; 
// #else 

  nat posA = rand.drawIntegerOpen(oldValues.size()) ; 
  nat posB = rand.drawIntegerOpen(oldValues.size()-1); 
  if(posA == posB)
    posB = oldValues.size() - 1 ; 

  // tout << "proposing " << posA << " and " << posB << std::endl ; 
  
  auto &rateA = oldValues.at(posA); 
  auto &rateB = oldValues.at(posB); 
  double sum = rateA + rateB; 

  // tout << "\n\noldVals="  << rateA << "," << rateB << std::endl; 
  
  double oldProb = rateA / (sum); 
  double newProb = rand.drawFromSlidingWindow(oldProb, parameter) ; 
  
  // tout << "new prob= " << newProb << std::endl; 

  if(newProb <= 0)
    newProb = - newProb; 
  else if(1 <= newProb )
    newProb = newProb-1; 

  // tout << "new prob after correct= " << newProb << std::endl; 
  
  assert(0. < newProb && newProb < 1. ); 

  rateA = newProb *   sum ; 
  rateB = (1-newProb) * sum ; 

  // tout << "newVals=" << rateA << ","  << rateB << std::endl; 

// #endif


  // check if contracts are met  
  if(minMaxIsRelative)
    {
      RateHelper::convertRelativeToLast(oldValues); 
      BoundsChecker::correctRevMat(oldValues); 
      RateHelper::convertToSum1(oldValues); 
    }
  else 
    {
      // sum may not be one 
      RateHelper::convertToSum1(oldValues) ; 
      correctAbsoluteRates(oldValues);
      RateHelper::convertToSum1(oldValues) ; 
    }

  if( fabs( std::accumulate(oldValues.begin(), oldValues.end(), 0.) - 1.0 ) > 1e-6 )
    {
      std::cerr << "Danger: while proposing values,  sum was " << std::accumulate(oldValues.begin(), oldValues.end(), 0.) << ". values: "   << oldValues  << std::endl; 
      assert(0); 
    }

  return oldValues; 
}



std::vector<double> SlidingProposer::proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
{
  if(oldValues.size() == 1 )
    return {proposeOneValue(oldValues[0], parameter,rand, hastings)}; 
  else
    return proposeRelativeMany(oldValues, parameter, rand, hastings);
}


