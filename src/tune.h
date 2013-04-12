#ifndef _SUCCESS_CTR
#define _SUCCESS_CTR


#include "axml.h"


/* #ifdef __cplusplus */
/* extern "C"{ */
/* #endif */


typedef struct 
{
  int lAcc; 			/// success within one tuning batch  
  int lRej; 
  
  int gAcc; 			/// overall acceptance / rejection 
  int gRej; 
} successCtr; 

void cntReject(successCtr *ctr); 
void cntAccept(successCtr *ctr); 
double getRatioOverall(successCtr *ctr); 
double getRatioLocal(successCtr *ctr); 
void resetCtr(successCtr *ctr); 

double tuneParameter(int batch, double accRatio, double parameter, boolean inverse ); 

/* #ifdef __cplusplus */
/* } */
/* #endif */
#endif
