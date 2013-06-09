#include "densities.h"
#include "axml.h"

#include <cassert>

using namespace std; 
#include <iostream>



static double betaFunctionLog(const vector<double> &alphas)
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



static double betaFunction(const vector<double> &alphas)
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



double gammaFunction(double alpha)
{
  
  double  lgam = lgamma(alpha);
  return signgam*exp(lgam); 
}




static void normalize( vector<double>  &vector, double normalizingConstant)
{
  double sum=0;
  for(auto v: vector)
    sum+=v; 
  
  for(auto &v : vector)
    v *= normalizingConstant/sum;   
}


double densityDirichletLog(vector<double> values, const vector<double> &alphas)
{
  double density=0;
  density -= betaFunctionLog(alphas);

  normalize(values, 1);

  for(int i=0; i<values.size(); i++)
    {      
      double val =  pow(values[i],alphas[i]-1); 
      density += log(val);
    }

  return density; 
}





double densityDirichlet(vector<double> values, const vector<double> &alphas)
{
  double density=1;

  density /= betaFunction(alphas);
  
  normalize(values, 1);
  
  for(int i=0; i< values.size(); i++)
    density *= pow( values[i],alphas[i]-1);
  
  return density; 
}


double exponentialDensity(double value, double lambda)
{
  return lambda * exp(- lambda * value); 
} 



