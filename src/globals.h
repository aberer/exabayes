#include "rng.h"

#ifdef _INCLUDE_DEFINITIONS

randKey_t rGlobalKey ; 
randCtr_t rGlobalCtr; 

/* more global variables =(  */
char configFileName[1024]; 
/* int seed;  */

char tree_file[1024]; 
char byteFileName[1024]; 

/* char** topologyFiles;  /\* a nexus-like output file (only sampled things)  *\/ */
/* char** outputParamFiles; /\*  outputs the parameters *\/ */
/* char infoFileName[1024];	   /\* outputs run specific information *\/ */


/* TODO  */
char binaryChainState[1024]; 

int numberOfStartingTrees = 0; 
int numberOfRuns = 0; 
int numberCoupledChains = 0; 


/* for crude performance measurements */
double timeIncrement = 0;  



#else 

/* global rng stuff  */
extern randKey_t rGlobalKey ; 
extern randCtr_t rGlobalCtr; 

extern int processID; 
extern char run_id[1024]; 
extern char configFileName[1024]; 
extern char workdir[1024]; 
extern char tree_file[1024]; 
extern char byteFileName[1024]; 


/* extern char** topologyFiles;  */
/* extern char** outputParamFiles;  */
extern char infoFileName[1024];

/* TODO */
extern char binaryChainState[1024]; 
extern int Thorough; 

extern int numberOfStartingTrees; 

extern int numberOfRuns; 
extern int numberCoupledChains; 


extern double timeIncrement;  
#endif
