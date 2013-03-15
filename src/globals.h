#include "rng.h"

#ifndef _GLOBALS_H
#define _GLOBALS_H


#include "config.h"

typedef struct globs
{
  int numberOfStartingTrees ; 
  int numberOfRuns; 
  int numberCoupledChains;

  /* RNG stuff */
  randKey_t rGlobalKey ; 
  randCtr_t rGlobalCtr; 

#if HAVE_PLL == 1 
  partitionList* partitions;   
#endif

  
  /* DEVEL not needed for prodctive runs  */
  state *allChains; 		/* careful with this! */
  int successFullSwitchesBatch; 

} globalAnalysisInfo; 
 

#endif


#ifdef _INCLUDE_DEFINITIONS

/* more global variables =(  */
char configFileName[1024]; 

char tree_file[1024]; 
char byteFileName[1024]; 

/* /\* TODO  *\/ */
/* char binaryChainState[1024];  */

globalAnalysisInfo gAInfo = {0,0,0}; 

/* for crude performance measurements */
double timeIncrement = 0;  

#else 

extern globalAnalysisInfo gAInfo; 

/* TODO  */
/* extern char binaryChainState[1024];  */

/* legacy */
extern char infoFileName[1024];
extern int Thorough; 
extern int processID; 
extern char run_id[1024]; 
extern char configFileName[1024]; 
extern char workdir[1024]; 
extern char tree_file[1024]; 
extern char byteFileName[1024]; 

extern double timeIncrement;  
#endif
