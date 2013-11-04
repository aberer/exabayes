
#include <cmath>
#include <cassert>
#include <iomanip>

#include "Branch.hpp"

#include "BoundsChecker.hpp"
#include "GlobalVariables.hpp"

// #include "TreeAln.hpp"

// here we initialize the max/min values for our various
// parameters. Static const is essentially like a global variable but
// saver. You access it with e.g., BoundsChecker::zmin, since the variable
// does not belong to an instance of the class, but the class in general. 

const double BoundsChecker::zMin = 1.0e-15 ; // 1e-16
const double BoundsChecker::zMax = (1.0 - 1.0e-6) ; // 1-1e-6
const double BoundsChecker::rateMin = 1e-7; 
const double BoundsChecker::rateMax = 1e6; 
const double BoundsChecker::alphaMin = 2e-2; 
const double BoundsChecker::alphaMax = 1e3; 
const double BoundsChecker::freqMin = 1e-3;


bool BoundsChecker::checkFrequencies( const std::vector<double> &freqs )  
{
  bool result = true; 
  double sum = 0; 
  for(auto &f : freqs)    
    {
      result &= freqMin <= f ; 
      sum += f; 
    }
  assert( fabs(sum  - 1.0 ) < 1e-6 ); 
  return result; 
}

bool BoundsChecker::checkBranch( const BranchLength &branch) 
{
  auto v = branch.getLength() ; 
  return BoundsChecker::zMin <= v && v <= BoundsChecker::zMax ; 
}

bool BoundsChecker::checkBranch( const BranchLengths &branch) 
{
  bool okay = true; 
  for(auto &v : branch.getLengths())
    okay &= BoundsChecker::zMin <= v && v <= BoundsChecker::zMax; 
  return okay; 
}

bool BoundsChecker::checkRevmat( const std::vector<double> &rates) 
{
  bool result = true; 
  for(auto &r : rates)
    {
      auto prob =  (r - BoundsChecker::rateMax < 1e-6 )  && ( BoundsChecker::rateMin <= r ); 
      if(not prob)
	tout << MAX_SCI_PRECISION << "attention: problem  with " << r << std::endl; 
      result &=  prob ; 
    }

  assert(*(rates.rbegin()) == 1.0); 
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


void BoundsChecker::correctBranch( BranchLengths &branch ) 
{
  auto lengths = branch.getLengths(); 

  for(auto &length : lengths )
    {
      if(length < BoundsChecker::zMin )
	length = BoundsChecker::zMin;
      if(BoundsChecker::zMax < length)
	length = BoundsChecker::zMax; 
    }

  branch.setLengths(lengths);   
} 


void BoundsChecker::correctBranch(BranchLength &branch)
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
  // because we have to convert the revmat from relative
  // representation to rate representation, this is more tricky.

  // note: rates are relative to last value here 

  // double sum = std::accumulate(rates.begin(), rates.end(), 0.); 

  for(auto &r : rates )
    {
      if( r <= BoundsChecker::rateMin)
	r = BoundsChecker::rateMin + std::numeric_limits<double>::epsilon();
      else if(BoundsChecker::rateMax <= r  ) 
	r = BoundsChecker::rateMax - std::numeric_limits<double>::epsilon();
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

