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

#include <thread>
#include <mutex>
#include "system/TeeStream.hpp"

class AdHocIntegrator; 
class TreeIntegrator; 

#define tout (*globals.teeOut) << SyncOut()

class TeeStream;

class GlobalVariables
{
public: 
  std::unique_ptr<TeeStream> teeOut; 
  std::unique_ptr<std::ofstream> logStream;
  std::mutex mtx;
};


#ifdef _INCLUDE_DEFINITIONS

std::ofstream nniOut; 
std::ofstream sprOut; 

GlobalVariables globals; 
std::chrono::system_clock::time_point timeIncrement;  
int debugPrint = 0; 
bool startIntegration = false; 

AdHocIntegrator* ahInt; 
TreeIntegrator* tInt; 

bool isYggdrasil; 

void (*exitFunction)(int code, bool waitForAll); 

std::thread::id _masterThread; 
volatile bool _threadsDie; 

#else 

extern std::thread::id _masterThread; 
extern volatile bool  _threadsDie; 

extern void (*exitFunction)(int code , bool waitForAll); 

extern bool isYggdrasil; 
extern int PLL_NUM_BRANCHES; 

extern std::ofstream nniOut; 
extern std::ofstream sprOut; 
extern bool startIntegration; 
extern AdHocIntegrator* ahInt; 
extern TreeIntegrator* tInt; 
extern GlobalVariables globals; 
extern std::chrono::system_clock::time_point timeIncrement;  
extern int debugPrint; 
#endif


#include "SyncOut.hpp"

#endif

