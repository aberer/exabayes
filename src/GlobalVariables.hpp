/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <string>
#include <chrono> 

#include "common.h"
#include "config.h"
#include "teestream.hpp"


// TODO this is still very bad style and we'll get rid of it as soon
// as I've the time to figure out how to do this 


// NOTICE it is not problem to change the numbers 


/* okay, so defining enums this way is rather save  */
#define NUM_PROPOSALS (19) //PROPOSALADD NUM_PROPOSALS NOTE Do not remove/modify  this line except for numerical value. The script addProposal.pl needs it as an identifier.
typedef enum
  {    
    ST_NNI= 0, 
    E_SPR = 1,
    E_TBR = 2,
    PARSIMONY_SPR= 3 , 
    GUIDED_SPR = 4,

    BRANCH_SLIDER = 5 ,
    TL_MULT = 6, 
    BRANCH_COLLAPSER = 7,
    NODE_SLIDER = 8,
    BRANCH_LENGTHS_MULTIPLIER = 9 , 
    UPDATE_SINGLE_BL_GUIDED = 10 , 

    REVMAT_SLIDER = 11 ,
    REVMAT_DIRICHLET = 12, 

    RATE_HET_SLIDER = 13, 
    RATE_HET_MULTI = 14, 
    
    FREQUENCY_SLIDER = 15, 
    FREQUENCY_DIRICHLET = 16, 
    
    AMINO_MODEL_JUMP = 17,

    BRANCH_GIBBS = 18

    //PROPOSALADD proposal_type NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
  } proposal_type;


#define tout ( *(globals.teeOut))
#define toutPure (*(globals.teeOut))

class TreeAln; 
class BipartitionHash; 

using namespace std; 


class GlobalVariables
{
public: 
  string logFile; 
  ofstream *logStream;   
  teestream* teeOut; 

#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTree; 
  bool verifyLnl;  		/* a hack around an ExaML problem. Just used for debugging */
#endif
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

#endif
