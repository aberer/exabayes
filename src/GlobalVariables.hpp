/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#ifndef _GLOBALS_H
#define _GLOBALS_H


#include "config.h"

class BipartitionHash; 
class SuccessCtr; 




class GlobalVariables
{
public: 
  analdef *adef; 

#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTree; 
  bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
#endif

  int myBatch ; 		/* if runs are executed in parallel: which runs should be done by this process? */
  int globalSize;
  int globalRank; 
};


#endif


#ifdef _INCLUDE_DEFINITIONS


bool isNewProposal[NUM_PROPOSALS] = 
  {
    true,     
    false, 
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    false,     
    true,     
    false,     
    false,     
    false,     
    false,     
    false,     
    true,
    true, 
    true 
  }; 


/* more global variables =(  */
char configFileName[1024]; 


GlobalVariables globals; 

/* for crude performance measurements */
double timeIncrement = 0;  

#else 
extern GlobalVariables globals; 

extern bool isNewPropossal[NUM_PROPOSALS]; 
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



