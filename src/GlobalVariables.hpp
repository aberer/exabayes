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

#include <fstream>
#include "common.h"
#include "config.h"
#include "teestream.hpp"

class AdHocIntegrator; 
class TreeIntegrator; 

#define tout (*globals.teeOut)


class GlobalVariables
{
public: 
  std::string logFile; 
  
  teestream* teeOut; 
  std::ofstream* logStream;
};

#endif

#ifdef _INCLUDE_DEFINITIONS

std::ofstream nniOut; 
std::ofstream sprOut; 

GlobalVariables globals; 
std::chrono::system_clock::time_point timeIncrement;  
int debugPrint = 0; 
bool startIntegration = false; 

AdHocIntegrator* ahInt; 
TreeIntegrator* tInt; 

#else 
extern std::ofstream nniOut; 
extern std::ofstream sprOut; 
extern bool startIntegration; 
extern AdHocIntegrator* ahInt; 
extern TreeIntegrator* tInt; 
extern GlobalVariables globals; 
extern int processID; 		// needed for raxml 
extern int processes; 
extern std::chrono::system_clock::time_point timeIncrement;  
extern int debugPrint; 
#endif
