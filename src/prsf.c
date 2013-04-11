
#include "axml.h" 
#include "bayes.h" 





/* double computePRSF(char *id, int numGen) */
/* { */
/*   printf("hello world\n");  */
/*   return 0;  */
/* } */













/**
   @brief the estimated potential scale rate reduction, as described by Yang.  
 */ 
double getPrsfForParameter(int paramNum, int numChain, int numGen, double ***matrix)
{
  double meanPerChain[numChain]; 
  
  for(int i = 0; i < numChain; ++i)
    {
      meanPerChain[i] = 0; 
      for(int j  = 0;  j < numGen; j++)
	{
	  meanPerChain[i] += matrix[i][paramNum][j]; 
	}
      meanPerChain[i] /= numGen; 
    }
  

  double overallMean = 0; 
  {
    for(int i = 0; i < numChain; ++i)
      overallMean += meanPerChain[i]; 
    overallMean /= numChain; 
  }

  double betweenChainVariance = 0; 
  {
    for(int i = 0; i < numChain; ++i)
      betweenChainVariance += pow(meanPerChain[i] - overallMean,2); 
    betweenChainVariance *= (numGen / (numChain-1)) ; 
  }


  double withinChainVariance = 0;    
  {
    double variancePerChain[numChain]; 

    for(int i = 0; i < numChain; ++i)
      {
	variancePerChain[i] = 0; 
	for(int j = 0; j < numGen; ++j)
	  variancePerChain[i] += pow(matrix[i][paramNum][j] -  meanPerChain[i],2);	  
	withinChainVariance += variancePerChain[i]; 
      }
    withinChainVariance /= ((numGen-1)* numChain); 
  }
  
  double tauSquare = ((numGen - 1)  / numGen ) * withinChainVariance  + (betweenChainVariance / numGen); 

  double prsf =  tauSquare / withinChainVariance; 
  
  return prsf; 
}



/**
   @brief parses the parameter output files 

   @param result -- a matrix of parameter values  mChains x parameter x nGen
   
 */ 
void parseFiles(int nGens, int mChains, int numParam, char *runId, double ****resultPtr, char ***paramNamesPtr)
{
  (*resultPtr) = exa_calloc(mChains, sizeof(double**)); 
  *paramNamesPtr = exa_calloc(numParam, sizeof(char*)); 
  for(int i = 0; i < numParam; ++i)
    (*paramNamesPtr)[i] = exa_calloc(1024, sizeof(char)); 

  /* scan for param names and determine which tokens to use */
  boolean* useToken = exa_calloc(2*numParam, sizeof(boolean)); 
  {
    char fullName[1024]; 
    sprintf(fullName, "ExaBayes_parameters.%s.%d", runId, 0); 
    FILE *fh = fopen(fullName, "r"); 
    size_t siz= 0; char *line = NULL; 
    getline(&line, &siz, fh); 	/* ignore fist line */
    getline(&line, &siz, fh) ;
    strtok(line, "\t"); 
    strtok(NULL, "\t"); 	/* forward by two  */
    char *pch = strtok(NULL, "\t");     
    int ctr = 0; 
    int nameCtr  = 0; 
    while(pch != NULL)
      {
	useToken[ctr] = strncmp(pch, "lnl", 3) !=  0; 	
	if(useToken[ctr])
	  {
	    strcpy((*paramNamesPtr)[nameCtr] , pch);
	    if((*paramNamesPtr)[nameCtr][strlen((*paramNamesPtr)[nameCtr])-1] == '\n')
	      (*paramNamesPtr)[nameCtr][strlen((*paramNamesPtr)[nameCtr])-1] = '\0';	      
	    nameCtr++; 
	  }
	ctr++; 
	pch = strtok(NULL, "\t");
      }    

    fclose(fh);
  } 

  for(int i = 0 ;i < mChains; ++i)
    {
      char fullName[1024]; 
      sprintf(fullName, "ExaBayes_parameters.%s.%d", runId, i); 
      FILE *fh = fopen(fullName, "r");       
      (*resultPtr)[i] = exa_calloc(numParam, sizeof(double*)); 
      for(int j = 0; j < numParam; ++j)
	(*resultPtr)[i][j] = exa_calloc(nGens, sizeof(double)); 

      size_t linesiz=0;
      char* line=NULL;

      /* ignore two lines */
      getline(&line, &linesiz, fh); 
      getline(&line, &linesiz, fh); 
      int sampNum = 0;
      while ( ( getline(&line, &linesiz, fh) > 0 )  && sampNum <nGens) 
	{
	  int ctr = 0, fillCtr = 0 ; 
	  char *pch = strtok(line, "\t"); 	  
	  pch = strtok(NULL, "\t"); 
	  pch  = strtok(NULL, "\t"); 
	  while(pch != NULL)
	  { 	    
	    if( useToken[ctr] )
	      (*resultPtr)[i][fillCtr++][sampNum] = atof(pch); 
	    ctr++; 	    
	    pch = strtok(NULL, "\t");
	  }
	  sampNum++; 
	} 
      fclose(fh); 
    }  

  
  exa_free(useToken);
}




void printMatrix(int numParam, int nGens, int mChains, char **paramNames, double ***result )
{
  for(int j = 0; j < numParam; ++j)
    {
      printf("%s\n", paramNames[j]);
      for(int k = 0; k < nGens; ++k)
  	{
  	  printf("%d\t", k);
  	  for(int i = 0; i < mChains; ++i)
  	    {
  	      printf("%.2f\t", result[i][j][k]);
  	    }
  	  printf("\n");
  	}
      printf("\n"); 
    }
}




static void freeMatrix(int mChains, int numParam, double ***result)
{
  for(int i = 0; i < mChains; ++i)
    {
      for(int j = 0; j < numParam; ++j)
	free(result[i][j]); 
      free(result[i]); 
    }
  free(result);
}



static int guessNumChain(char *runid)
{
  boolean exists = TRUE; 
  int ctr = 0; 

  do 
    {
      char fullName[1024]; 
      sprintf(fullName, "ExaBayes_parameters.%s.%d", runid,ctr); 
      /* printf("trying to open %s\n", fullName);  */
      FILE *fh = fopen(fullName, "r"); 
      if(fh == NULL)
	{	  
	  exists = FALSE; 
	}
      else 
	{
	  fclose(fh); 
	  ctr++;       
	}
    } while(exists); 
  
  assert(ctr  > 0 ); 
  return ctr; 
}



static int guessNumParam(char *runid)
{
  int ctr = 0; 
  
  char fullName[1024]; 
  sprintf(fullName, "ExaBayes_parameters.%s.%d", runid, ctr);   
  
  FILE *fh = fopen(fullName, "r"); 
  assert(fh != NULL); 
  size_t siz = 0; 
  char *line = NULL;  
  getline(&line, &siz, fh);   
  assert(siz); 
  getline(&line, &siz, fh); 	/* second line  */
  assert(siz); 

  char *pch = strtok(line, "\t"); 
  pch = strtok(NULL, "\t"); 
  pch = strtok(NULL, "\t"); 
  
  while(pch != NULL)
    {
      if(strncmp(pch, "lnl",3) != 0 )
	ctr++; 
      pch = strtok(NULL, "\t"); 
    }
  
  fclose(fh); 
  return ctr; 
}


/**
   @brief guesses the numebr of generations that have been finished by ALL chains 
   
   @return number of samples that have been completed by EVERY chain
   
   Notice, we have to know the number of parameters first in order to
   check, if the lines have been written appropriately
 */
int guessNumGen( int numChain, int numParam, char *runid)
{
  int maxCompleted = 0; 
  
  for(int i = 0; i < numChain; ++i)
    {
      char fullName[1024]; 
      sprintf(fullName, "ExaBayes_parameters.%s.%d", runid, i); 
      FILE *fh = fopen(fullName,"r"); 
      assert(fh); 
      char *line = NULL; size_t siz = 0; 
      getline(&line,&siz,fh ); 
      getline(&line,&siz,fh ); 

      int numTokens = 0; 
      {
	char *pch = strtok(line, "\t"); 
	while(pch != NULL)
	  {	    
	    pch = strtok(NULL, "\t") ;
	    numTokens++; 
	  }
      }

      int goodLines = 0; 
      while( getline(&line, &siz, fh) > 0)
	{
	  int paramsSeen = 0; 
	  char *pch = strtok(line, "\t");
	  while(pch != NULL)
	    {
	      paramsSeen++; 
	      pch = strtok(NULL, "\t"); 
	    }
	  if(numTokens == paramsSeen)
	    ++goodLines;
	} 

      maxCompleted = (i == 0) ? goodLines : MIN(goodLines, maxCompleted); 
    }

  return maxCompleted; 
}




/**
   @brief reads all parameter files and prints the PRSF for each
   parameter

   @notice NOT thread-safe at the moment! 
 */ 

double printPRSF(char *runId )
{  
  double ***matrix = NULL; 
  char **names = NULL; 

  int numChain =  guessNumChain(runId), 
    numParam =  guessNumParam(runId),
    numGen = guessNumGen(numChain, numParam, runId); 

  parseFiles(numGen,numChain, numParam, "testRun",&matrix, &names);

  printf("PRSF by component:\n"); 
  int window = 10; 
  for(int iter = 0 ; iter < numParam / window; ++iter)
    {
      for(int i = iter * window; i < (iter+1) * window && i < numParam ; ++i)
	printf("%s\t", names[i]);     
      printf("\n"); 
      for(int i = iter * window ;  i < (iter+1) * window && i < numParam; ++i)
	{
	  double localPrsf = getPrsfForParameter(i, numChain, numGen, matrix); 
	  printf("%.3f\t\t",  localPrsf);
	}	    
      printf("\n\n"); 
    }

  freeMatrix(numChain, numParam, matrix);
}





#ifdef _STANDALONE
int main(int argc, char **argv)
{
  if(argc < 2 || 4 < argc)
    {
      printf("USAGE: %s runId\nMust be called from folder where the respective files are.\n ", argv[0]); 
      exit(-1); 
    }

  char runId[1024]; 
  strcpy(runId, argv[1]); 

  
  double ***matrix = NULL; 
  char **names = NULL; 

  int numChain =  guessNumChain(runId), 
    numParam =  guessNumParam(runId),
    numGen = guessNumGen(numChain, numParam, runId); 
  
  printf("found %d chains with %d parameters and %d good lines \n", numChain, numParam,numGen );

  parseFiles(numGen,numChain, numParam, "testRun",&matrix, &names);

  printf("printing parsed samples by chain:\n");
  printMatrix(numParam, numGen, numChain, names, matrix);

  printf("PRSF for component:\n"); 
  for(int i = 0; i < numParam;++i)
    {
      double localPrsf = getPrsfForParameter(i, numChain, numGen, matrix); 
      printf("%s\t\t%.3f\n", names[i], localPrsf);
    }

  freeMatrix(numChain, numParam, matrix); 
  return 0;     
}
#endif





