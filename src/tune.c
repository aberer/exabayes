

#include <math.h>
#include "tune.h" 

#include "axml.h"
#include "bayes.h"  
#include "globals.h"








/**
   @brief log tuning for a parameter 
   
   @return tuned parameter    
 */
double tuneParameter(int batch, double accRatio, double parameter, boolean inverse )
{  
  double delta = 1.0 / sqrt(batch);
  delta = 0.01 < delta ? 0.01 : delta;

  double logTuning = log(parameter);
  
  if(inverse)
    logTuning += (accRatio > TARGET_RATIO)  ? -delta : +delta ;
  else 
    logTuning += (accRatio > TARGET_RATIO)  ? +delta : -delta ;
  
  double newTuning = exp(logTuning);

  /* TODO min+max tuning?  */
  /* if (newTuning > minTuning && newTuning < maxTuning) */

  return newTuning; 
}





/**
   @brief resets the local counter.
*/ 
void resetCtr(successCtr *ctr)
{
  ctr->lAcc = 0; 
  ctr->lRej = 0; 
}



/**
    @brief Gets the local success rate  
 */
double getRatioLocal(successCtr *ctr)
{
  return (double)ctr->lAcc / ((double)(ctr->lAcc + ctr->lRej)) ; 
} 


/**
    @brief Gets the overall success rate  
 */
double getRatioOverall(successCtr *ctr)
{
  return (double)ctr->gAcc / ((double)(ctr->gAcc + ctr->gRej)) ; 
}


void cntAccept(successCtr *ctr)
{
  ctr->lAcc++; 
  ctr->gAcc++; 
}


void cntReject(successCtr *ctr)
{
  ctr->lRej++; 
  ctr->gRej++; 
}
