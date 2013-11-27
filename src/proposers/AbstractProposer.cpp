#include <algorithm>
#include "AbstractProposer.hpp"

AbstractProposer::AbstractProposer(bool tune, bool tuneup, double minVal, double maxVal)
  : _tune{tune}
  , _tuneup{tuneup}
  , _minVal{minVal}
  , _maxVal{maxVal}
{
  
} 


void AbstractProposer::correctAbsoluteRates(std::vector<double> &values) const 
{
  nat prevFixed = -1; 
  nat nowFixed = 0; 
  
  auto fixedHigh = std::vector<bool>(values.size(), false);
  auto fixedLow = std::vector<bool>(values.size(), false);

  nat iter = 0; 
  bool valuesNotOkay = true; 
  while(prevFixed != nowFixed && valuesNotOkay)
    {
      double normer = 0; 
      nat ctr = 0; 
      for(auto &v : values)
	{
	  assert(v > 0); 	    
	  if( v <= _minVal)
	    {
	      // tout << "================ value too small!" << std::endl; 
	      v = _minVal; 
	      fixedLow.at(ctr) = true; 
	    }
	  else if(_maxVal <= v)
	    {
	      // tout << "================ value too high!" << std::endl; 
	      v = _maxVal; 
	      fixedHigh.at(ctr) = true; 
	    }
	  else 
	    normer += v; 
	  ++ctr ; 
	}

      nat numHigh = std::count_if(begin(fixedHigh), end(fixedHigh), [](bool elem) {return elem; }); 
      nat numLow = std::count_if(begin(fixedLow), end(fixedLow), [](bool elem){return elem; }); 

      normer = ( 1. - (_minVal * numLow  +  _maxVal * numHigh))    / normer  ;

      
      valuesNotOkay = false; 
      ctr = 0; 
      for(auto &v : values)
	{
	  if(not  ( fixedHigh.at(ctr)  ||  fixedLow.at(ctr))  )
	    {
	      v *= normer; 

	      // meh: could simplify everything alltogether again 
	      valuesNotOkay |= ( v <= _minVal || _maxVal <= v) ; 
	    }
	  ++ctr; 
	}

      prevFixed = nowFixed; 
      nowFixed = numHigh + numLow; 
      ++iter; 
    }

} 

