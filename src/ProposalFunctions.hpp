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
    updateHastings(hastings, densityDirichlet(oldValues, newValues) / densityDirichlet(newValues,oldValues) , "dirichlet"); 

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
    updateHastings(hastings, multiplier, "multiplier"); 
    oldValues[position] = oldValues[position] * multiplier; 
    return oldValues;
  }

  static bool tune; 
  static bool tuneup; 
}; 



class SlidingProposal
{
public: static vector<double> getNewValues(vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {
    int position = rand.drawRandInt(oldValues.size());
    oldValues[position] = rand.drawFromSlidingWindow(oldValues[position], parameter); 

    // this is only ONE way to do it! 
    double sum = 0; 
    for(auto v : oldValues)
      sum += v ; 

    for(auto &v : oldValues)
      v /= sum ; 

    return oldValues; 
  }  

  static bool tune; 
  static bool tuneup; 
}; 


#endif
