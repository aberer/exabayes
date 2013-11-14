#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "Density.hpp"
#include "axml.h"

// bracen copy from mrb
static double logGamma (double alp)
{
  double 
    x = alp, f = 0.0, z;
	
  if (x < 7) 
    {
      f = 1.0;
      z = x-1.0;
      while (++z < 7.0)  
	f *= z;
      x = z;   
      f = -log(f);
    }
  z = 1.0 / (x*x);
  return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
	   (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
	    0.083333333333333)/x);  

}


static double betaFunctionLog(const std::vector<double> &alphas)
{
  double beta=0;
  double sum=0;
  for(auto v : alphas)
    {
      beta+=lgamma(v);       
      sum+=v ;    
    }
  beta -= lgamma(sum);

  return beta;
}


static double gammaFunction(double alpha)
{
  
  double  lgam = lgamma(alpha);
  return signgam*exp(lgam); 
}


static double betaFunction(const std::vector<double> &alphas)
{
  double beta=1.0;
  double sum=0;

  for(auto v : alphas)
    {
      beta *= gammaFunction(v); 
      sum += v; 
    }

  beta /= gammaFunction(sum);

  return beta;
}

namespace Density
{
  double lnDirichlet(std::vector<double> values, const std::vector<double> &alphas)
  {
    double density=0;
    density -= betaFunctionLog(alphas);
    
    // double sum = 0; 
    if( fabs( std::accumulate(begin(values), end(values), 0. )    - 1. )  > 1e-6)
      {
	double sum = std::accumulate(begin(values), end(values), 0. ); 
	std::cout << "error: expected sum of 1, got " << sum << std::endl; 
	assert(0); 
      }      

    for(nat i=0; i<values.size(); i++)
      density += (alphas[i] - 1 ) * log(values[i]); 

    return density; 
  }

  double dirichlet(std::vector<double> values, const std::vector<double> &alphas)
  {
    assert(0); 			// deprecated 
    double density=1;  
    density /= betaFunction(alphas);  
    assert(fabs(std::accumulate(values.begin(), values.end(), 0.) -  1.0 )  < 1e-6); 

    for(nat i=0; i< values.size(); i++)
      density *= pow( values[i],alphas[i]-1);
  
    return density; 
  }


  double lnExponential(double value, double lambda)
  {
    return  log(lambda) - lambda * value; 
  }

  double exponential(double value, double lambda)
  {
    assert(0); 			// deprecated 
    return lambda * exp(- lambda * value); 
  } 

  double lnGamma (double x, double alpha, double beta )
  {    
    double result =  (alpha-1.0) * log(x) + alpha * log(beta) - x * beta - logGamma(alpha);  

    // std::cout << "Gamma("<< x << ";" <<  alpha << "," << beta << ") = " << result << std::endl; 
    return result; 
  }
}  
