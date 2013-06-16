#include "axml.h" 
#include "tune.h"

/**
   @brief log tuning for a parameter 
   
   @return tuned parameter    
 */
double tuneParameter(int batch, double accRatio, double parameter, bool inverse )
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
  
  double minTuning = 1e-8,
    maxTuning = 1e5; 
  if (minTuning <  newTuning && newTuning < maxTuning)
    return  newTuning; 
  else 
    return parameter; 
}
