#ifndef _PROPOSAL_FUNCTION_H
#define _PROPOSAL_FUNCTION_H

#include <vector>

#include "Randomness.hpp"
#include "densities.h"
#include "Chain.hpp"
#include "AbstractProposal.hpp"

using namespace std; 

class DirichletProposal
{				
public: static vector<double> getNewValues(vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {
    double sum = 0; 
    for(auto &v : oldValues)
      sum += v; 
    assert(fabs(sum - 1.0 ) < 1e-6); 

    vector<double> newValues; 
    rand.drawDirichletExpected(newValues, oldValues, parameter * oldValues.size() );
    AbstractProposal::updateHastings(hastings, densityDirichlet(oldValues, newValues) / densityDirichlet(newValues,oldValues) , "dirichlet"); 

    sum = 0; 
    for(auto &v : newValues)
      sum += v; 
    assert(fabs(sum - 1.0)< 1e-6); 

    return newValues; 
  }
  
  static bool tune; 
  static bool tuneup; 
}; 


class ExponentialProposal
{
public: static vector<double> getNewValues(vector<double> oldValue, double parameter, Randomness &rand, double &hastings)
  {
    // TODO @kassian: how to modify the hastings? 
    // assert(0); 
    return vector<double>();

  }

  static bool tune; 
  static bool tuneup; 
  
};


class BiunifProposal
{  
  // TODO incorporate the parametr 
public: static vector<double> getNewValues(vector<double> oldValue, double parameter, Randomness &rand, double &hastings)
  {
    // TODO @ kassian: how to modify the hastings?  
    assert(0); 
    return vector<double>();
  }
  
  static bool tune; 
  static bool tuneup; 
};  



class MultiplierProposal
{
public: static vector<double> getNewValues(vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {    
    int position = rand.drawRandInt(oldValues.size()); 
    double multiplier =  rand.drawMultiplier( parameter); 
    AbstractProposal::updateHastings(hastings, multiplier, "multiplier"); 
    oldValues[position] = oldValues[position] * multiplier; 
    return oldValues;
  }

  static bool tune; 
  static bool tuneup; 
}; 



#define ALT_SLIDER  


class SlidingProposal
{
public: static vector<double> getNewValues(vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {
    if(oldValues.size() == 1 )
      {
	double newVal = rand.drawFromSlidingWindow(oldValues[0], parameter);
	if(newVal < 0 )
	  newVal = - newVal ; 
	oldValues[0] = newVal; 
      }
    else
      {
#ifdef ALT_SLIDER
	int posA = rand.drawRandInt(oldValues.size()); 
	double oldVal = oldValues[posA] ; 
	double newVal = rand.drawFromSlidingWindow(oldVal, parameter); 


#if TODO
	// use the min values here 
	assert(0); 
#endif

	if(newVal < 0 )
	  newVal = -newVal; 
	if(newVal > 1) 
	  newVal =  newVal - 1 ; // TODO cannot be good 

	oldValues[posA] = newVal; 


	double sum = 0; 
	for(auto v : oldValues)
	  {
	    sum  += v; 
	    assert(v > 0); 	    
	  }
	
	for(auto &v : oldValues)
	  v /= sum; 
	  
#else 
	int posA = rand.drawRandInt(oldValues.size()),
	  posB = rand.drawRandInt(oldValues.size()-1); 
	if(posB == posA)
	  posB = oldValues.size()-1; 

	double both =  (oldValues[posA] + oldValues[posB]); 
	double oldProp = oldValues[posA] / both ; 
	double newProp = rand.drawFromSlidingWindow(oldProp, parameter);
	if(newProp < 0 )
	  newProp = - newProp; 
	if(newProp > 1 )
	  newProp = newProp - 1 ; 

	oldValues[posA] = newProp * both; 
	oldValues[posB] = (1-newProp) * (both); 

	double sum = 0; 
	for(auto v : oldValues)
	  {
	    sum += v ; 
	    assert(v > 0) ; 
	  }
	
	assert(fabs( sum - 1.0 ) < 1e-6);       
#endif
	
      }

    return oldValues; 
  }  

  static bool tune; 
  static bool tuneup; 
}; 


#endif
