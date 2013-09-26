#include <algorithm>
#include "AbstractProposer.hpp"


void AbstractProposer::correctAbsoluteRates(std::vector<double> &values) const 
{
  nat prevFixed = -1; 
  nat nowFixed = 0; 
  
  auto fixedHigh = std::vector<bool>(values.size(), false);
  auto fixedLow = std::vector<bool>(values.size(), false);
  
  nat iter = 0; 
  while(prevFixed != nowFixed)
    {
      double normer = 0; 
      nat ctr = 0; 
      for(auto &v : values)
	{
	  assert(v > 0); 	    
	  if( v <= minVal)
	    {
	      v = minVal; 
	      fixedLow[ctr] = true; 
	    }
	  else if(maxVal <= v)
	    {
	      v = maxVal; 
	      fixedHigh[ctr] = true; 
	    }
	  else 
	    normer += v; 
	  ++ctr ; 
	}

      nat numHigh = std::count_if(fixedHigh.begin(), fixedHigh.end(), [](bool elem) {return elem; }); 
      nat numLow = std::count_if(fixedLow.begin(), fixedLow.end(), [](bool elem){return elem; }); 

      normer = ( 1 - (minVal * numLow  +  maxVal * numHigh))    / normer  ;

      ctr = 0; 
      for(auto &v : values)
	{
	  if(not  fixedHigh[ctr]  && not  fixedLow[ctr] )
	    v *= normer; 
	  ++ctr; 
	}

      prevFixed = nowFixed; 
      nowFixed = numHigh + numLow; 
      ++iter; 
    }
} 

