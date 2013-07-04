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


bool BoundsChecker::checkFrequencies( const std::vector<double> &freqs )  
{
  return true; 
}

bool BoundsChecker::checkBranch( const Branch &branch ) 
{
  return BoundsChecker::zMin <= branch.getLength() && branch.getLength()  <= BoundsChecker::zMax ; 
}

bool BoundsChecker::checkRevmat( const std::vector<double> &rates) 
{
  bool result = true; 
  for(auto &r : rates)
    result &=  (  r  <=  BoundsChecker::rateMax)  && ( BoundsChecker::rateMin <= r ) ; 
  return result; 
}
 
bool BoundsChecker::checkAlpha( double alpha) 
{
  return BoundsChecker::alphaMin <= alpha    && alpha <= BoundsChecker::alphaMax; 
} 


void BoundsChecker::correctAlpha(double &alpha) 
{
  if( alpha < BoundsChecker::alphaMin)
    alpha = BoundsChecker::alphaMin; 
  if(  BoundsChecker::alphaMax < alpha )
    alpha = BoundsChecker::alphaMax; 
}


void BoundsChecker::correctBranch( Branch &branch ) 
{
  double length = branch.getLength(); 
  if(length < BoundsChecker::zMin )
    length = BoundsChecker::zMin;
  if(BoundsChecker::zMax < length)
    length = BoundsChecker::zMax; 
  branch.setLength(length); 
}
 
void BoundsChecker::correctRevMat( std::vector<double> &rates)
{
  for(auto &r : rates )
    {
      if( r < BoundsChecker::rateMin)
	r = BoundsChecker::rateMin; 
      else if(BoundsChecker::rateMax < r  )
	r = BoundsChecker::rateMax; 
    }  
}
 
void BoundsChecker::correctFrequencies( std::vector<double> &frequencies)
{
  // determine number of freqs that are not okay 
  int numberOkay = 0; 
  double sum = 0;  
  for(auto &f : frequencies)
    {
      if(f < BoundsChecker::freqMin) 
	f = BoundsChecker::freqMin; 
      else 
	{
	  ++numberOkay; 
	  sum += f ; 
	}
    }

  // renormalize again 
  for(auto &f : frequencies)
    if(f != BoundsChecker::freqMin) 
	f /= sum; 
  
  // check to be sure that they add up to 1 
  sum = 0 ;
  for(auto &f : frequencies)
    sum += f; 
  assert(abs(sum - 1.0 ) < 1e-6); 
} 
