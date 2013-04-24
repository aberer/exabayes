#ifndef _MAIN_COMMON_H 
#define _MAIN_COMMON_H 





#if HAVE_PLL != 0
void initializeTree(tree *tr, partitionList *partitions, analdef *adef); 
#else  
void initializeTree(tree *tr, analdef *adef); 
#endif

/* provides  */
void printVersionInfo(); 
void printREADME(); 

/* int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg); */ 
/* void get_args(int argc, char *argv[], analdef *adef, tree *tr);  */
void parseCommandLine(int argc, char *argv[], analdef *adef); 
void analyzeRunId(char id[128]); 
/* void initAdef(analdef *adef);  */
void ignoreExceptionsDenormFloat(); 
int filexists(char *filename); 
void finalizeFiles(); 

#endif
