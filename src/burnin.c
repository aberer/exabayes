#include "axml.h" 
#include "bayes.h"


typedef struct 
{
  double value; 
  int origPos; 
} rankHelper; 


static int sortByPostProb(const void *a, const void *b)
{
  return ((rankHelper*)a)->value - ((rankHelper*)b)->value; 
}




static double getSpearmanInWindow(double *values, int length)
{
  rankHelper *ranks = exa_calloc(length, sizeof(rankHelper)); 
  for(int i = 0; i < length; ++i)
    {
      ranks[i].value = values[i]; 
      ranks[i].origPos = i; 
    }

  qsort(ranks, length,sizeof(rankHelper), sortByPostProb); 
  
  double sum = 0; 
  for(int i = 0; i < length; ++i)
    sum += pow(ranks[i].origPos - i,2) ; 
  sum *= 6 ; 
  sum /= length * ((length * length) - 1  ); 
  
  double result = 1 - sum; 
  exa_free(ranks); 
  
  return result; 
}

  




/**
   @brief recognizes burnin by correlation of posterior probability
   values

   @return the first generation where we are sure about having reached
   stationarity
 
   @param pPValues -- the posterior probabilities by sampled generation 
   @param interval -- sampling interval
 */
int recognizeBurninBySpearmanCorrelation(double *pPValues, int length, int windowSize)
{
  for(int i = 0; i < length / windowSize; ++i)    
    {
      double* curStart = pPValues + (i * windowSize); 
      double cor = getSpearmanInWindow(curStart, windowSize); 
      int start = (i * (length / windowSize)); 
      int end = ((i+1) * (length / windowSize)); 
      printf("%d-%d: %.2f\n", start, end, cor); 
    }
  
  return 0 ; 
}
