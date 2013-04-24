#include "axml.h"
#include "bayes.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

// #include "main-common.h"
#include "prsfComputer.h"

#if HAVE_PLL != 0 
#include "globalVariables.h" 
#endif

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
