/**
    @file globals.h

    @brief Global variables for ExaBayes.
    
    Notice that there are already a bunch of global variables from axml-variants.  
*/ 


#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <string>
#include "config.h"
#include "teestream.hpp"


// TODO this is still very bad style and we'll get rid of it as soon
// as I've the time to figure out how to do this 

#define NUM_PROP_CATS 5 
typedef enum _cats {
  TOPOLOGY = 1  , 
  BRANCH_LENGTHS = 2, 
  FREQUENCIES = 3,
  SUBSTITUTION_RATES = 4  ,
  RATE_HETEROGENEITY = 5
} category_t; 



/* okay, so defining enums this way is rather save  */
#define NUM_PROPOSALS (24) //PROPOSALADD NUM_PROPOSALS NOTE Do not remove/modify  this line except for numerical value. The script addProposal.pl needs it as an identifier.
typedef enum
  {
    TL_MULT = 0, 
    UPDATE_MODEL = 1 ,
    UPDATE_GAMMA = 2 ,
    UPDATE_GAMMA_EXP=3,
    UPDATE_SINGLE_BL = 4,
    UPDATE_SINGLE_BL_EXP = 5 ,
    UPDATE_SINGLE_BL_BIUNIF = 6,
    UPDATE_MODEL_BIUNIF = 7,
    UPDATE_MODEL_SINGLE_BIUNIF = 8,
    UPDATE_MODEL_ALL_BIUNIF = 9,
    UPDATE_MODEL_PERM_BIUNIF = 10,
    UPDATE_FREQUENCIES_BIUNIF = 11,
    BRANCH_LENGTHS_MULTIPLIER =12,
    FREQUENCY_SLIDER = 13,
    GUIDED_SPR = 14,
    ST_NNI= 15, 
    GAMMA_MULTI = 16,
    NODE_SLIDER = 17,
    UPDATE_FREQUENCIES_DIRICHLET = 18,
    UPDATE_MODEL_DIRICHLET = 19,
    UPDATE_SINGLE_BL_GUIDED = 20,
    E_SPR = 21,
    E_TBR = 22,
    PARSIMONY_SPR= 23 
    //PROPOSALADD proposal_type NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
  } proposal_type;








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


// i really hate to this ... 
// category_t categoriesOfProposals = 
//   {
//     TOPOLOGY, 
//     SUBSTITUTION_RATES, 
//     RATE_HETEROGENEITY, 
//     RATE_HETEROGENEITY, 

    
//   } ; 



GlobalVariables globals; 

/* for crude performance measurements */
double timeIncrement = 0;  

#else 
extern GlobalVariables globals; 
extern int processID; 		// needed for raxml 
extern double timeIncrement;  

#endif






