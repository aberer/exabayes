/**************************************************************************************/
/* NOTICE: contrary to its name, this file cannot be shared yet across pll and examl  */
/**************************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include "axml.h"
#include "main-common.h"
#include "globals.h" 

int modelExists(char *model, tree *tr);
int filexists(char *filename);  


int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
}


void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    bad_opt    =FALSE,
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",          
    *optarg,
    model[1024] = ""
  ;

  double 
    likelihoodEpsilon;
  
  int     
    optind = 1,        
    c,
    nameSet = 0,
    treeSet = 0,   
    modelSet = 0, 
    byteFileSet = 0;


  /*********** tr inits **************/ 
 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;

  tr->manyPartitions = FALSE;

  tr->categories             = 25;

  tr->grouped = FALSE;
  tr->constrained = FALSE;

  tr->gapyness               = 0.0; 
  tr->saveBestTrees          = 0;

  tr->useMedian = FALSE;

  seed = -1 ; 

  strcpy(configFileName, "\0");

  /********* tr inits end*************/




  while(!bad_opt && ((c = mygetopt(argc,argv,"e:c:f:i:m:t:w:n:s:vhMQap:", &optind, &optarg))!=-1))
    {
    switch(c)
      {    
      case 'a':
	tr->useMedian = TRUE;
	break;
      case 'Q':
	tr->manyPartitions = TRUE;   	
	break;
      case 's': 
	strcpy(byteFileName, optarg);	 	
	byteFileSet = TRUE;
	break;      
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;                                 
      case 'e':
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;
	break;          
      case 'v':
	/* TODO */
	/* printVersionInfo(); */
	errorExit(0);
      
      case 'h':
	/* TODO */
	/* printREADME(); */
	errorExit(0); 
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;
      case 'n':
        strcpy(run_id,optarg);
	/* TODO */
	/* analyzeRunId(run_id); */
	nameSet = 1;
        break;
      case 'w':
        strcpy(resultDir, optarg);
	resultDirSet = TRUE;
        break;
      case 't':
	strcpy(tree_file, optarg);       
	treeSet = 1;       
	break;     
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, tr) == 0)
	  {
	    if(processID == 0)
	      {
		printf("Rate heterogeneity Model %s does not exist\n\n", model);               
		printf("For per site rates (called CAT in previous versions) use: PSR\n");	
		printf("For GAMMA use: GAMMA\n");		
	      }
	    errorExit(-1);
	  }
	else
	  modelSet = 1;
	break;     
      case 'c': 
	strcpy(configFileName, optarg); 
	break; 
      case 'p':
	seed = atoi(optarg); 
	break; 
      default:
	errorExit(-1);
      }
    }

  
  if(strlen(configFileName) == 0 ||  ! filexists(configFileName) )
    {
      if(processID == 0)
	printf("\nPlease provide a config file via -C <file>. A testing file is available under  \n"); 
      errorExit(-1); 
    }

  
  if(seed == -1 )
    {
      if(processID == 0)
	printf("\nPlease provide a proper seed for the initialization of the random number generator. \n "); 
      errorExit(-1); 
    }

  if( ! filexists( configFileName))
    {
      if(processID == 0 )
  	printf("\nPlease provide a minimal config file via -C <file>.\n");
      errorExit(-1);
    }


  if(!byteFileSet)
    {
      if(processID == 0)
	printf("\nError, you must specify a binary format data file with the \"-s\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    {
      if(processID == 0)
	printf("\nError, you must specify a model of rate heterogeneity with the \"-m\" option\n");
      errorExit(-1);
    }

  if(!nameSet)
    {
      if(processID == 0)
	printf("\nError: please specify a name for this run with -n\n");
      errorExit(-1);
    }

  if(!treeSet && !adef->useCheckpoint)
    {
      if(processID == 0)
	{
	  printf("\nError: please either specify a starting tree for this run with -t\n");
	  printf("or re-start the run from a checkpoint with -R\n");
	}
      
      errorExit(-1);
    }
  
   {

    const 
      char *separator = "/";

    if(resultDirSet)
      {
	char 
	  dir[1024] = "";
	

	if(resultDir[0] != separator[0])
	  strcat(dir, separator);
	
	strcat(dir, resultDir);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	strcpy(workdir, dir);
      }
    else
      {
	char 
	  dir[1024] = "",
	  *result = getcwd(dir, sizeof(dir));
	
	assert(result != (char*)NULL);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	
	strcpy(workdir, dir);		
      }
   }

  return;
}
