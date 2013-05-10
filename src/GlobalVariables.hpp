/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#ifndef _GLOBALS_H
#define _GLOBALS_H


#include "config.h"
#include "proposalType.h"

class BipartitionHash; 
class SuccessCtr; 


class GlobalVariables
{
public: 
#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTree; 
  bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
#endif
};


#endif


#ifdef _INCLUDE_DEFINITIONS

bool isNewProposal[NUM_PROPOSALS] = 
  {
    true,     
    true, 
    true,     
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
    true,     
    false,     
    true,     
    true,     
    false,     
    true,     
    true,     
    false,     
    true,
    true, 
    true 
  }; 


/* more global variables =(  */
// char configFileName[1024]; 


GlobalVariables globals; 

/* for crude performance measurements */
double timeIncrement = 0;  

#else 
extern GlobalVariables globals; 

extern bool isNewPropossal[NUM_PROPOSALS]; // HACK  needed 
extern int processID; 		// needed for raxml 
extern double timeIncrement;  

// extern char run_id[1024]; 
// extern char configFileName[1024]; 
// extern char workdir[1024]; 
// extern char tree_file[1024]; 
// extern char byteFileName[1024]; 
extern char infoFileName[1024];

#endif
