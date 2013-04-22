/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#include "rng.h"

#ifndef _GLOBALS_H
#define _GLOBALS_H


#include "config.h"

class BipartitionHash; 
class SuccessCtr; 


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

  state *allChains; 		/* careful with this! */
  int samplingFrequency; 
  /* hashtable *bvHash;  */
  BipartitionHash* bipHash; 

  int diagFreq;   
  int numGen;			/// just relevent, if we have exactly 1 run 


  
  double* temperature;   /// one temperature for each run 
  /* successCtr **swapInfo;   /// indicates how ofter swaps  between chains succeeded -> only the upper half w/o diag is filled! */
  
  /* todo replace with a chain swapper? */
  SuccessCtr **swapInfo; 
  
  analdef *adef; 

#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTree; 
#endif
  int tuneFreq; 
  int printFreq; 
  double asdsfIgnoreFreq; 
  double asdsfConvergence; 
  double heatFactor; 
  int swapInterval; 
  bool tuneHeat;
  int burninGen ; 
  double burninProportion; 
  
  bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
} globalAnalysisInfo; 

 

#endif


#ifdef _INCLUDE_DEFINITIONS

/* more global variables =(  */
char configFileName[1024]; 


globalAnalysisInfo gAInfo = 
  {0,
   1,
   1,
   {{0,0}},
   {{0,0}},
   NULL, 
   0,
   NULL,
   0,
   10000 ,
   NULL , 
   NULL,  
   NULL, 
#ifdef DEBUG_LNL_VERIFY
   NULL, 
#endif   
   100, 
   500, 
   0.1, 
   0.005, 
   0.1 ,
   1,
   true, 
   5000, 
   0.0,
   true   
  }; 


/* for crude performance measurements */
double timeIncrement = 0;  

#else 

extern globalAnalysisInfo gAInfo; 



extern int Thorough; 
extern int processID; 
extern char run_id[1024]; 
extern char configFileName[1024]; 
extern char workdir[1024]; 
extern char tree_file[1024]; 
extern char byteFileName[1024]; 
extern double timeIncrement;  
extern char infoFileName[1024];

#endif



