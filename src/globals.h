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

  int samplingFrequency; 
  /* BipartitionHash* bipHash;  */

  int diagFreq;   
  int numGen;			/// just relevent, if we have exactly 1 run 

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


globalAnalysisInfo gAInfo; 

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



