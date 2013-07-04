#include "BoundsChecker.hpp"


// here we initialize the max/min values for our various
// parameters. Static const is essentially like a global variable but
// saver. You access it with e.g., BoundsChecker::zmin, since the variable
// does not belong to an instance of the class, but the class in general. 

const double BoundsChecker::zMin = 1.0E-15 ; 
const double BoundsChecker::zMax = (1.0 - 1.0E-6) ; 
const double BoundsChecker::rateMin = 0.0000001; 
const double BoundsChecker::rateMax = 1000000.0; 
const double BoundsChecker::alphaMin = 0.02; 
const double BoundsChecker::alphaMax = 1000.0; 
const double BoundsChecker::freqMin = 0.001; 


bool checkFrequencies( const std::vector<double> &freqs )  
{
  return true; 
}

bool checkBranch( const Branch &branch ) 
{
  return BoundsChecker::zMin < branch.getLength() && branch.getLength()  < BoundsChecker::zMax ; 
}

bool checkRevmat( const std::vector<double> &rates) 
{
  bool result = true; 
  for(auto &r : rates)
    result &=  (  r  <  BoundsChecker::rateMax)  && ( BoundsChecker::rateMin < r ) ; 
  return result; 
}
 
bool checkAlpha( double alpha) 
{
  return BoundsChecker::alphaMin < alpha    && alpha < BoundsChecker::alphaMax; 
} 


