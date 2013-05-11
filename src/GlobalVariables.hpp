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


#define tout (*(globals.teeOut))


class TreeAln; 
class BipartitionHash; 

using namespace std; 


class GlobalVariables
{
public: 

#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTree; 
  bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
#endif
  string logFile; 
  ofstream *logStream;   
  teestream* teeOut; 
};




#endif


#ifdef _INCLUDE_DEFINITIONS

GlobalVariables globals; 

/* for crude performance measurements */
double timeIncrement = 0;  

#else 
extern GlobalVariables globals; 
extern int processID; 		// needed for raxml 
extern double timeIncrement;  

#endif
