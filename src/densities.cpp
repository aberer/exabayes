#include "densities.h"
#include "axml.h"

#include <cassert>



double betaFunction(double *alpha, int length)
{
  double beta=1.0;
  double sum=0;
  for(int i=0; i<length;i++)
    {
      beta*=gammaFunction(alpha[i]); 
      sum+=alpha[i];    
    }
  beta=beta/(gammaFunction(sum));
 
  return beta;
}





double gammaFunction(double alpha)
{
  
  double  lgam = lgamma(alpha);
  return signgam*exp(lgam); 
}




static void normalize(double * vector, int length, double normalizingConstant)
{
  double sum=0;
  for(int i=0; i<length; i++)
    sum+=vector[i];
  
  for(int i=0; i<length; i++)
    vector[i]=vector[i]*normalizingConstant/sum;
  
}



/** 
    @brief just wraps the method for usage with vectors 
 */ 
double densityDirichletWrapper(vector<double> values, vector<double> alphas)
{
  assert(values.size() == alphas.size());
  double tmpVal[values.size()] , 
    tmpAlpha[alphas.size()]; 
  for(nat i = 0; i < values.size(); ++i)
    {
      tmpVal[i] = values[i]; 
      tmpAlpha[i] = alphas[i]; 
    }

  return densityDirichlet(tmpVal, tmpAlpha, values.size());
}



double densityDirichlet(double *values, double *alphas, int length)
{
  double density=1;
  double normValues[length];
  
  for(int i=0; i<length; i++)
    normValues[i]=values[i];
  
  density=density/betaFunction(alphas, length);
  
  normalize(normValues, length, 1);
  
  for(int i=0; i<length; i++)
    density=density*pow(normValues[i],alphas[i]-1);
  
  return density; 
}


double exponentialDensity(double value, double lambda)
{
  return lambda * exp(- lambda * value); 
} 



