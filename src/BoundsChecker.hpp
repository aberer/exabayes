#ifndef _BOUNDS_CHECKER
#define _BOUNDS_CHECKER

#include <vector>

#include "Branch.hpp"


class BoundsChecker
{
public: 
  static const double zMin, zMax, 
    rateMax, rateMin, 
    alphaMin, alphaMax,
    freqMin; 

  static bool checkFrequencies( const std::vector<double> &freqs )  ; 
  static bool checkBranch( const Branch &branch ) ; 
  static bool checkRevmat( const std::vector<double> &rates) ; 
  static bool checkAlpha( double alpha) ; 
  

  static void correctAlpha(double &alpha) ;   
  static void correctBranch( Branch &branch ) ; 
  static void correctRevMat( std::vector<double> &rates); 
  static void correctFrequencies( std::vector<double> &frequencies); 

} ; 


#endif
