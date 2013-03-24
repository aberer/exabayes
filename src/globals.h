/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#include "rng.h"

#ifndef _GLOBALS_H
#define _GLOBALS_H


#include "config.h"


/* Those globals are bad.  Having them as a singleton is slightly
   better. 
   
   NOTICE Mind the initialization! 

 */
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
  state *allChains; 		/* careful with this! */
  int samplingFrequency; 
  hashtable *bvHash; 
  int diagFreq;   
  int numGen;			/// just relevent, if we have exactly 1 run 

  accRejCtr **swapInfo;   /// indicates how ofter swaps  between chains succeeded -> this is a half-matrix for each independent run

} globalAnalysisInfo; 

 

#endif


#ifdef _INCLUDE_DEFINITIONS

/* more global variables =(  */
char configFileName[1024]; 

char tree_file[1024]; 
char byteFileName[1024]; 

/* /\* TODO  *\/ */
/* char binaryChainState[1024];  */

globalAnalysisInfo gAInfo = 
  {0,
   0,
   0,
   {{0,0}},
   {{0,0}},
   #if HAVE_PLL == 1 
   NULL,
   #endif
   NULL, 
   0,
   NULL,
   0,
   10000   ,
   NULL,
  }; 


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
