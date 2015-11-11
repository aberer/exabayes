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

class TreeAln; 
class BipartitionHash; 

class GlobalVariables
{
public: 
  std::string logFile; 
  
  teestream* teeOut; 
  std::ofstream* logStream;

};


#endif

#ifdef _INCLUDE_DEFINITIONS


GlobalVariables globals; 
std::chrono::system_clock::time_point timeIncrement;  

#else 
extern GlobalVariables globals; 
extern int processID; 		// needed for raxml 
extern std::chrono::system_clock::time_point timeIncrement;  

#endif
