

#ifdef _INCLUDE_DEFINITIONS

/* more global variables =(  */
char configFileName[1024]; 
int seed; 

#else 

extern int seed; 
extern char byteFileName[1024]; 
extern char binaryCheckpointInputName[1024]; 
extern int processID; 
extern char run_id[1024]; 
extern char tree_file[1024]; 
extern char configFileName[1024]; 
extern char workdir[1024]; 

#endif
