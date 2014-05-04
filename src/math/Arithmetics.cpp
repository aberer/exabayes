#include "Arithmetics.hpp"
#include "common.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <cassert>


namespace Arithmetics
{
  double getCoefficientOfVariation(const std::vector<double> &data) 
  {
    return sqrt(getVariance(data)) / getMean(data); 
  } 


  double getPercentile(double percentile, std::vector<double> data)
  {
    assert(percentile < 1.);
    assert(data.size() > 0 ); 

    if(data.size() == 1)
      return data.at(0); 

    sort(begin(data), end(data), std::less<double>()); 

    nat idx = nat(double(data.size()) * percentile); 

    if( fabs(double(idx) / percentile - data.size()) > 0   ) // dirty
      ++idx;

    if(data.size() < idx + 1 )	// meh 
      idx = data.size() -1 ; 
    
    auto result = data.at(idx); 
    // return ; 
    assert(not std::isinf(result) );
    return result; 
  }

  double getMean(const std::vector<double> &data)
  {
    double result = 0; 
    for(auto &v : data)
      result += v; 

    result /= data.size(); 
    
    return result; 
  }

  double getVariance(const std::vector<double> &data)
  {
    if(data.size() == 1)
      return 0 ; 

    double result = 0; 
    auto mean = Arithmetics::getMean(data);
    
    for(auto &v : data)
      result += pow(mean - v,2 ); 
    result /= (data.size() -1);  
    
    return result; 
  }
  
  // this implementation is inspired by:
  // https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDAQFjAA&url=http%3A%2F%2Fwww.people.fas.harvard.edu%2F~plam%2Fteaching%2Fmethods%2Fconvergence%2Fconvergence_print.pdf&ei=J_0iUrWWE8STtQb8uoHoDg&usg=AFQjCNE6hl7PsS1PX6vU3nQjWQNNaEBIiQ&sig2=yq_0iw4pw7svqMS054dDbQ&cad=rja
  double PRSF(const std::vector<std::vector<double> > data)
  {
    if(data.size() == 1)
      return 1; 

    double withinChainVariance = 0;   
    for(auto &chainData : data)
      withinChainVariance += getVariance(chainData); 
    withinChainVariance /= data.size(); 

    auto means = std::vector<double>{}; 
    for(auto &chainData: data ) 
      means.push_back(Arithmetics::getMean(chainData)); 
    double meanOfMeans = Arithmetics::getMean(means); 
    double betweenChainVariance = 0; 

    // slight deviation: the chains may not have the same number of
    // draws (i.e., branch lengths). So for the between chain variance
    // (bcv), we do not multiply the overall bcv with n (as in the
    // formula), but each component gets multiplied with the number of
    // draws in the chain.
    // if we have the same number of draws, this is correct anyway 
    // TODO sounds reasonable, but must be verified

    auto numberOfDraws = std::vector<double>{}; 
    for(auto &chainData : data)
      numberOfDraws.push_back(double(chainData.size())); 
    auto meanNumberOfDraws = Arithmetics::getMean(numberOfDraws);
    // the above sounds less good   

    for(nat i = 0; i < means.size() ;++i )
      betweenChainVariance += data[i].size() * pow(means[i] - meanOfMeans, 2); 
    betweenChainVariance /= (data.size()-1) ; 

    double estimatedVariance = 0; // of the stationary chain
    estimatedVariance = (1. - 1. / meanNumberOfDraws) * withinChainVariance + betweenChainVariance / double(meanNumberOfDraws) ;  
    
    return sqrt(estimatedVariance / withinChainVariance); 
  }
  
  
  double getPearsonCorrelationCoefficient(const std::vector<double> &sampleA, const std::vector<double> &sampleB)
  {
    assert(sampleA.size() == sampleB.size()); 
    
    auto meanA = getMean(sampleA) ; 
    auto meanB = getMean(sampleB); 
    
    auto varA = getVariance(sampleA); 
    auto varB = getVariance(sampleB); 

    double result = 0; 
    for(nat i = 0; i < sampleA.size() ;++i)
      result += (sampleA[i] - meanA) * (sampleB[i] - meanB) ; 
    result /= (sqrt(varA) * sqrt(varB)); 
    
    return result; 
  }

  
  // implementation according to   

  // makes some sense 
  // http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introbayes_sect008.htm

  double getAutoCorrelation(const std::vector<double> &data, nat lag)
  {
    auto mean = getMean(data); 

    double denominator = 0.; 
    for(auto v : data )
      denominator += (v - mean) * (v - mean); 

    double result = 0.f; 
    for(int i = 0; i < int(data.size() - lag); ++i)
      result += ( data[i] - mean ) * (data[i+lag] - mean) ; 
    
    result /= denominator; 
    return result; 
  }


  double getEffectiveSamplingSize(const std::vector<double> &data)
  {
    if(data.size() == 0)
      return 0; 
    if(data.size() == 1)
      return 1; 

    double gammaStat[2000];
    double mean = getMean(data) ; 
    nat maxLag = data.size()-1; 

    if(data.size()  - 1 > 2000)
      maxLag = 2000; 
    double varStat = 0; 
    
    for (nat lag = 0; lag < maxLag; lag++)
      {
	gammaStat[lag]=0;
	for (nat j = 0; j < data.size() - lag; j++) 
	  {
	    double del1 = data[j] - mean;
	    double del2 = data[j + lag] - mean;
	    gammaStat[lag] += (del1 * del2);
	  }

	gammaStat[lag] /= ((double) (data.size() - lag));

	if (lag == 0) 
	  varStat = gammaStat[0];
	else if (lag % 2 == 0) 
	  {
	    if (gammaStat[lag - 1] + gammaStat[lag] > 0) 
	      varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
	    else
	      maxLag = lag;
	  }
      }
    return (gammaStat[0] * data.size()) / varStat;
  }
}
