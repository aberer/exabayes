/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <string>
#include <chrono> 
#include <memory>

#include "common.h"
#include "config.h"
#include "teestream.hpp"

#define tout (*globals.teeOut)

// #define tout ( *(globals.teeOut))
// #define toutPure (*(globals.teeOut))

class TreeAln; 
class BipartitionHash; 

using namespace std; 


class GlobalVariables
{
public: 
  string logFile; 
  
  teestream* teeOut; 
  ofstream* logStream;
  


// #ifdef DEBUG_LNL_VERIFY
//   TreeAln *debugTree; 
//   bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
// #endif
};


#endif

#ifdef _INCLUDE_DEFINITIONS


GlobalVariables globals; 
/* for crude performance measurements */
// double timeIncrement = 0;  
chrono::system_clock::time_point timeIncrement;  

#else 
extern GlobalVariables globals; 
extern int processID; 		// needed for raxml 
extern chrono::system_clock::time_point timeIncrement;  
// extern unique_ptr<teestream> teeOut; 
// extern unique_ptr<ofstream> logStream ;


#endif
