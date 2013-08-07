#include <cmath>

#include <iomanip>
#include "BoundsChecker.hpp"
#include "GlobalVariables.hpp"

#include "TreeAln.hpp"

// here we initialize the max/min values for our various
// parameters. Static const is essentially like a global variable but
// saver. You access it with e.g., BoundsChecker::zmin, since the variable
// does not belong to an instance of the class, but the class in general. 

const double BoundsChecker::zMin = 1.0e-15 ; // 1e-16
const double BoundsChecker::zMax = (1.0 - 1.0e-6) ; // 1-1e-6
const double BoundsChecker::rateMin = 1e-7; 
const double BoundsChecker::rateMax = 1e6; 
const double BoundsChecker::alphaMin = 0.02; 
const double BoundsChecker::alphaMax = 1000.0; 
const double BoundsChecker::freqMin = 0.001;


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

bool BoundsChecker::checkBranch( const Branch &branch) 
{
  bool okay = true; 
  for(auto &v : branch.getAllLengths())
    {
      // TODO what a hack!  
      if(v != TreeAln::problematicBL)
	okay &= BoundsChecker::zMin <= v && v <= BoundsChecker::zMax; 
    }
  return okay; 
}

bool BoundsChecker::checkRevmat( const std::vector<double> &rates) 
{
  bool result = true; 
  for(auto &r : rates)
    result &=  (  r  <=  BoundsChecker::rateMax)  && ( BoundsChecker::rateMin <= r ) ; 
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


void BoundsChecker::correctBranch(Branch &branch, const AbstractParameter* param)
{
  double length = branch.getLength(param); 
  if(length < BoundsChecker::zMin )
    length = BoundsChecker::zMin;
  if(BoundsChecker::zMax < length)
    length = BoundsChecker::zMax; 
  branch.setLength(length, param);   
}


void BoundsChecker::correctBranch( Branch &branch, const std::vector<AbstractParameter*> &params ) 
{
  for(auto &param : params)
    BoundsChecker::correctBranch(branch, param); 
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
