#include <algorithm>
#include "AbstractProposer.hpp"


void AbstractProposer::correctAbsoluteRates(std::vector<double> &values) const 
{
  // tout << MAX_SCI_PRECISION << "to be normalized " << values << std::endl; 


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
	      // tout << "value too small!" << std::endl; 
	      v = minVal; 
	      fixedLow.at(ctr) = true; 
	    }
	  else if(maxVal <= v)
	    {
	      // tout << "value too high!" << std::endl; 
	      v = maxVal; 
	      fixedHigh.at(ctr) = true; 
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
	  if(not  ( fixedHigh.at(ctr)  ||  fixedLow.at(ctr))  )
	    v *= normer; 
	  ++ctr; 
	}

      prevFixed = nowFixed; 
      nowFixed = numHigh + numLow; 
      ++iter; 
    }

  // tout << "iterations necessary: " <<  iter << "\tminval="<< minVal << std::endl; 
} 

