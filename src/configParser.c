#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "axml.h"
#include "proposals.h"
#include "globals.h"


extern void initDefaultValues(state *theState); 


int parseLine(char* line, int length, state *theState)
{
  char
    key[1024], 
    value[1024],
    *tok; 

  tok = strtok(line, "="); 
  if(tok == NULL)
    return 1; 
  sscanf(tok,  "%s", key); 

  tok = strtok(NULL, "="); 
  if(tok == NULL) 
    return 1; 
  sscanf(tok, "%s", value); 

  if( ! strcmp(key, "initSPRWeight"))
    theState->proposalWeights[SPR] = atof(value); 
  else if( ! strcmp(key, "initModelWeight"))
    theState->proposalWeights[UPDATE_MODEL] = atof(value); 
  else if( ! strcmp(key, "initGammaWeight"))
    theState->proposalWeights[UPDATE_GAMMA] = atof(value); 
  else if( ! strcmp(key, "initSingleBranchWeight"))
    theState->proposalWeights[UPDATE_SINGLE_BL] = atof(value);   
  else if( ! strcmp(key, "initSingleBranchExpWeight"))
    theState->proposalWeights[UPDATE_SINGLE_BL_EXP] = atof(value);   
  else if( ! strcmp(key, "numGen"))
    theState->numGen = atoi(value);   
  else if( ! strcmp(key, "initPenaltyFactor"))
    theState->penaltyFactor = atof(value); 
  else 
    {
      printf("could not parse option %s in line %s\n", tok, line); 
      assert(0); 
    }

  return 0; 
}


void parseConfig(state *theState)
{
  ssize_t read = 0; 
  size_t length = 0; 
  int ctr = 1; 
  char *line = NULL; 
  
  initDefaultValues(theState);

  FILE *fh = myfopen(configFileName, "r"); 
  
  while( (read = getline(&line, &length, fh) )   != -1 )
    {
      if(parseLine(line, length,  theState) > 0)
	{
	  printf("Error while parsing config file. Could not process line number %d\n", ctr);
	  exit(-1);
	}
      ctr++;
    }
}

