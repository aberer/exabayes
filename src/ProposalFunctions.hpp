
#ifndef _PROPOSAL_FUNCTION_H
#define _PROPOSAL_FUNCTION_H

#include "Randomness.hpp"
#include "densities.h"
#include <vector>

using namespace std; 

class DirichletProposal
{				
public: static vector<double> getNewValues(vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {
    vector<double> newValues; 
    double tmp[oldValues.size()]; 
    for(nat i = 0; i < oldValues.size(); ++i)
      tmp[i] = oldValues[i]; 
    double tmpNew[oldValues.size()];
    rand.drawDirichletExpected(tmpNew, tmp, parameter, (int)oldValues.size());
    for(nat i = 0; i < oldValues.size(); ++i)
      newValues.push_back( tmpNew[i]); 
    
    // cout << "old :" ; 
    // for(auto v : oldValues)
    //   cout << v << "," ; 
    // cout << endl; 

    // cout << "new :" ; 
    // for(auto v : newValues)
    //   cout << v << ",";
    // cout << endl; 

    // cout << oldValues << endl; 
    // cout << newValues << endl; 


    hastings *= densityDirichlet(tmp, tmpNew, oldValues.size()) / densityDirichlet(tmpNew,tmp, oldValues.size());
    
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
    hastings *= multiplier; 
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
    // no hastings moification needed
    return oldValues; 
  }  

  static bool tune; 
  static bool tuneup; 
}; 


#endif
