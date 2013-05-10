/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <string>
#include "config.h"
#include "proposalType.h"
#include "teestream.hpp"


class TreeAln; 
class BipartitionHash; 
class SuccessCtr; 

using namespace std; 


class GlobalVariables
{
public: 

#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTree; 
  bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
#endif
  string logFile; 
  teestream* tout; 
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


GlobalVariables globals; 

/* for crude performance measurements */
double timeIncrement = 0;  

#else 
extern GlobalVariables globals; 

extern bool isNewPropossal[NUM_PROPOSALS]; // HACK  needed 
extern int processID; 		// needed for raxml 
extern double timeIncrement;  

#endif
