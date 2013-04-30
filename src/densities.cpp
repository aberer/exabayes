#include "densities.h"

#include "axml.h"



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
