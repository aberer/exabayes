#ifndef _PROPOSAL_FUNCTION_H
#define _PROPOSAL_FUNCTION_H

#include <vector>
#include <algorithm>

#include "Randomness.hpp"
#include "densities.h"
#include "Chain.hpp"
#include "AbstractProposal.hpp"

//////////////
// ABSTRACT //
//////////////
class AbstractProposer
{
public:   
  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings) = 0; 

  bool isTune() const {return tune; } 
  bool isTuneup() const {return tuneup; }

  // TODO would be cool 
  // void setTunedParameter(double param ) { tunedParameter = param; }

protected: 
  bool tune; 
  bool tuneup; 
  double minVal; 
  double maxVal; 
}; 



///////////////
// DIRICHLET //
///////////////
class DirichletProposal : public AbstractProposer
{				
public: 
  DirichletProposal( double minVal, double maxVal) 
  {
    this->minVal = minVal; 
    this->maxVal = maxVal; 
    tune = true ; 
    tuneup = false; 
  } 

  
  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {
    assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 

    std::vector<double> newValues; 
    
    std::vector<double> scaledOld = oldValues; 
    for_each(scaledOld.begin(), scaledOld.end(), [&](double &d) { return d *= parameter * oldValues.size() ; }) ; 
    newValues = rand.drawRandDirichlet(scaledOld); 
    
    // newValues = rand.drawDirichletExpected( oldValues, parameter * oldValues.size() );

    // correct for problematic values 
    int numHigh = 0,
      numLow = 0; 
    double normer = 0; 
    for(auto &v : newValues) 
      {
	if(v < minVal)
	  {
	    v = minVal; 
	    ++numLow ; 
	  }
	else if(maxVal < v)
	  {
	    v = maxVal; 
	    ++numHigh; 
	  }
	else 
	  normer += v ; 
      }

    normer =   (1. -( numHigh *maxVal + numLow * minVal))  / normer; 
    
    for(auto &v : newValues)
      {
	if( v != minVal && v != maxVal)
	  v *= normer; 
      }

    AbstractProposal::updateHastings(hastings, densityDirichlet(oldValues, newValues) / densityDirichlet(newValues,oldValues) , "dirichlet"); 
    
    assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0. ) - 1.0 ) < 1e-6); // sum of rates equals 1 ? 
    return newValues; 
  }

}; 


////////////////
// MULTIPLIER //
////////////////
class MultiplierProposal : public AbstractProposer
{
public: 

  MultiplierProposal(double minVal, double maxVal)
  {
    this->minVal = minVal; 
    this->maxVal = maxVal; 
    tune = true; 
    tuneup = true; 
  }

  virtual std::vector<double> proposeValues(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {    
    double newVal = 0, position = 0, multiplier = 0; 

    position = rand.drawRandInt(oldValues.size()); // dont like 
    multiplier =  rand.drawMultiplier( parameter);     
    newVal = oldValues[position] * multiplier; 

    // TODO allowed? 
    if(newVal < minVal)
      newVal = minVal; 
    else if(maxVal < newVal)
      newVal = maxVal; 

    AbstractProposal::updateHastings(hastings, multiplier, "multiplier"); 
    oldValues[position] = newVal; 
    return oldValues;
  }
}; 


#define ALT_SLIDER  

//////////////
// SLIDER   //
//////////////
class SlidingProposal : public AbstractProposer
{
public: 
  SlidingProposal(double minVal, double maxVal)
  {
    this->minVal = minVal; 
    this->maxVal = maxVal; 
    tune = true; 
    tuneup = true; 
  }


  double proposeOneValue(double oldVal, double parameter, Randomness &rand, double &hastings)
  {
    double newVal = rand.drawFromSlidingWindow(oldVal, parameter);
    if(newVal < 0 )
      newVal = - newVal ; 
    assert(0); 			// TODO  
    return newVal; 
  }


  // important: this is under the assumption, that stuff sums up to one 
  std::vector<double> proposeRelativeMany(std::vector<double> oldValues, double parameter, Randomness &rand, double &hastings)
  {
#ifdef ALT_SLIDER
    int posA = rand.drawRandInt(oldValues.size()); 

    double oldVal = oldValues[posA] ; 
    double newVal = rand.drawFromSlidingWindow(oldVal, parameter); 

    if(newVal < 0 )
      newVal = -newVal; 
    if(newVal > 1) 
      newVal =  newVal - 1 ; // TODO cannot be good 

    oldValues[posA] = newVal; 

    double normer = 0; 
    nat numHigh = 0,
      numLow = 0; 
    for(auto &v : oldValues)
      {
	assert(v > 0); 	    
	if( v < minVal)
	  {
	    v = minVal; 
	    ++numLow; 
	  }
	else if(maxVal < v)
	  {
	    v = maxVal; 
	    ++numHigh;
	  }
	else 
	  normer += v; 
      }
    
    normer = (1 - (minVal * numLow  +  maxVal * numHigh)) / normer ;
    // if(oldValues.size()  == 6 )
    //   std::cout << normer << std::endl; 

    for(auto &v : oldValues)
      if(v != minVal && v != maxVal)
	v *= normer; 

    assert(fabs(std::accumulate(oldValues.begin(), oldValues.end(), 0.)  - 1.0 ) < 1e-6 ); 
    
    // if(oldValues.size() == 4)
    //   {
    // 	for_each(oldValues.begin(), oldValues.end(), [] (double &d ){tout << std::setprecision(4)  << d << "," ; }) ; 
    // 	tout << std::endl; 
    //   }
#else 

    // NOT WORKING
    assert(0); 

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
    return oldValues; 
  }


  virtual std::vector<double> proposeValues(std::vector<double> oldValues, 
					    double parameter, Randomness &rand, double &hastings)
  {
    if(oldValues.size() == 1 )
      return {proposeOneValue(oldValues[0], parameter,rand, hastings)}; 
    else
      return proposeRelativeMany(oldValues, parameter, rand, hastings);
  }  
}; 


#endif





// disabling those for now. Just more stuff to maintain. 
// we can easily reactivate them 
#if 0 
class ExponentialProposal : public AbstractProposer
{
public: 
  ExponentialProposal()
  {
    tune = false; 
    tuneup = false; 
  }

  virtual std::vector<double> proposeValues(std::vector<double> oldValue, double parameter, Randomness &rand, double &hastings)
  {
    // TODO @kassian: how to modify the hastings? 
    // assert(0); 
    return std::vector<double>();

  }

};


class BiunifProposal : public AbstractProposer
{  
  // TODO incorporate the parametr 
public: 
  BiunifProposal()
  {
    tune = true; 
    tuneup = true; 
  }

  virtual std::vector<double> proposeValues(std::vector<double> oldValue, double parameter, Randomness &rand, double &hastings)
  {
    // TODO @ kassian: how to modify the hastings?  
    assert(0); 
    return std::vector<double>();
  }

} ;  

#endif
