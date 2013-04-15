#include "axml.h"
#include "bayes.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "main-common.h"

#include "proposals.h"
#include "output.h"
/* turn on, when in release mode */
/* #define PRODUCTIVE */

#include "adapters.h"


#include "globalVariables.h" 

int processID = 0; 
int seed; 


void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile); 

static void printBoth(FILE *f, const char* format, ... )
{
  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);
}



static void printModelAndProgramInfo(tree *tr, partitionList *pr, analdef *adef, int argc, char *argv[])
{
  FILE *infoFile = myfopen(infoFileName, "ab");
  char modelType[128];


  if(tr->useMedian)
    strcpy(modelType, "GAMMA with Median");
  else
    strcpy(modelType, "GAMMA");


  printVersionInfo();

  if(!adef->compressPatterns)
    printBoth(infoFile, "\nAlignment has %d columns\n\n",  tr->originalCrunchedLength);
  else
    printBoth(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->originalCrunchedLength);

  printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");

  if(adef->perGeneBranchLengths)
    printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", pr->numberOfPartitions);
  else
    printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", pr->numberOfPartitions);

  for( int model = 0; model < pr->numberOfPartitions; model++)
  {
    printBoth(infoFile, "Partition: %d\n", model);
    printBoth(infoFile, "Alignment Patterns: %d\n", pr->partitionData[model]->upper - pr->partitionData[model]->lower);
    printBoth(infoFile, "Name: %s\n", pr->partitionData[model]->partitionName);

    switch(pr->partitionData[model]->dataType)
    {
      case DNA_DATA:
        printBoth(infoFile, "DataType: DNA\n");
        printBoth(infoFile, "Substitution Matrix: GTR\n");
        break;
      case AA_DATA:
        assert(pr->partitionData[model]->protModels >= 0 && pr->partitionData[model]->protModels < NUM_PROT_MODELS);
        printBoth(infoFile, "DataType: AA\n");
        printBoth(infoFile, "Substitution Matrix: %s\n", protModels[pr->partitionData[model]->protModels]);
        printBoth(infoFile, "%s Base Frequencies:\n", (pr->partitionData[model]->protFreqs == 1)?"Empirical":"Fixed");
        break;
      case BINARY_DATA:
        printBoth(infoFile, "DataType: BINARY/MORPHOLOGICAL\n");
        printBoth(infoFile, "Substitution Matrix: Uncorrected\n");
        break;
      case SECONDARY_DATA:
        printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");
        printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
        break;
      case SECONDARY_DATA_6:
        printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");
        printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
        break;
      case SECONDARY_DATA_7:
        printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");
        printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
        break;
      case GENERIC_32:
        printBoth(infoFile, "DataType: Multi-State with %d distinct states in use (maximum 32)\n",pr->partitionData[model]->states);
        switch(tr->multiStateModel)
        {
          case ORDERED_MULTI_STATE:
            printBoth(infoFile, "Substitution Matrix: Ordered Likelihood\n");
            break;
          case MK_MULTI_STATE:
            printBoth(infoFile, "Substitution Matrix: MK model\n");
            break;
          case GTR_MULTI_STATE:
            printBoth(infoFile, "Substitution Matrix: GTR\n");
            break;
          default:
            assert(0);
        }
        break;
      case GENERIC_64:
        printBoth(infoFile, "DataType: Codon\n");
        break;
      default:
        assert(0);
    }
    printBoth(infoFile, "\n\n\n");
  }

  printBoth(infoFile, "\n");

  printBoth(infoFile, "%s was called as follows:\n\n", PROGRAM_NAME);
  for( int i = 0; i < argc; i++)
    printBoth(infoFile,"%s ", argv[i]);
  printBoth(infoFile,"\n\n\n");

  fclose(infoFile);
}


void readByteFile(tree *tr, analdef *adef, partitionList *partitions, double ***empiricalFrequencies)
{
  unsigned char *y;

  FILE 
    *byteFile = myfopen(byteFileName, "rb");	 

  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
  myBinFread(&(partitions->numberOfPartitions),  sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);

  partitions->perGeneBranchLengths = adef->perGeneBranchLengths;

  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed RAxML-specific */

  tr->aliaswgt                   = (int *)exa_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
  myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       

  tr->rateCategory    = (int *)    exa_malloc((size_t)tr->originalCrunchedLength * sizeof(int));	  

  tr->patrat          = (double*)  exa_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->patratStored    = (double*)  exa_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
  tr->lhs             = (double*)  exa_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 

  *empiricalFrequencies = (double **)exa_malloc(sizeof(double *) * partitions->numberOfPartitions);

  y = (unsigned char *)exa_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));
  tr->yVector = (unsigned char **)exa_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

  for( nat i = 1; i <= (size_t)tr->mxtips; i++)
    tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength]; 

  setupTree(tr, FALSE, partitions);

  for(int i = 0; i < partitions->numberOfPartitions; i++)
    partitions->partitionData[i]->executeModel = TRUE;

  /* data structures for convergence criterion need to be initialized after! setupTree */

  for( nat i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int len;
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr->nameList[i] = (char*)exa_malloc(sizeof(char) * (size_t)len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      /*printf("%s \n", tr->nameList[i]);*/
    }  

  for( nat i = 1; i <= (size_t)tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  for(int model = 0; model < partitions->numberOfPartitions; model++)
    {
      int 
	len;

      pInfo 
	*p = partitions->partitionData[model];

      myBinFread(&(p->states),             sizeof(int), 1, byteFile);
      myBinFread(&(p->maxTipStates),       sizeof(int), 1, byteFile);
      myBinFread(&(p->lower),              sizeof(int), 1, byteFile);
      myBinFread(&(p->upper),              sizeof(int), 1, byteFile);
      myBinFread(&(p->width),              sizeof(int), 1, byteFile);
      myBinFread(&(p->dataType),           sizeof(int), 1, byteFile);
      myBinFread(&(p->protModels),         sizeof(int), 1, byteFile);
      myBinFread(&(p->autoProtModels),     sizeof(int), 1, byteFile);
      myBinFread(&(p->protFreqs),          sizeof(int), 1, byteFile);
      myBinFread(&(p->nonGTR),             sizeof(boolean), 1, byteFile);
      myBinFread(&(p->numberOfCategories), sizeof(int), 1, byteFile); 

      /* later on if adding secondary structure data

	 int    *symmetryVector;
	 int    *frequencyGrouping;
      */

      myBinFread(&len, sizeof(int), 1, byteFile);
      p->partitionName = (char*)exa_malloc(sizeof(char) * (size_t)len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);

      (*empiricalFrequencies)[model] = (double *)exa_malloc(sizeof(double) * (size_t)partitions->partitionData[model]->states);
      myBinFread((*empiricalFrequencies)[model], sizeof(double), partitions->partitionData[model]->states, byteFile);
    }

  myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

  fclose(byteFile);
}




void initializeTree(tree *tr, partitionList *partitions, analdef *adef)
{
  partitions->partitionData = (pInfo**)exa_malloc(NUM_BRANCHES*sizeof(pInfo*));
  double **empiricalFrequencies;
  readByteFile(tr, adef, partitions, &empiricalFrequencies); 
  initializePartitionsSequential(tr, partitions);
  initModel(tr, empiricalFrequencies, partitions);
}



void exa_main(tree *tr, 
#if HAVE_PLL == 1 
	      partitionList *partitions,
#endif
	      analdef *adef); 

int main (int argc, char *argv[])
{ 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* parallel stuff will not be easy to set up  */
  assert(0); 
#endif

  ignoreExceptionsDenormFloat(); 
  tree  *tr = (tree*)exa_calloc(1,sizeof(tree));
  analdef *adef = (analdef*)exa_calloc(1,sizeof(analdef));

  masterTime = gettime();         
  initAdef(adef);
  get_args(argc, argv, adef, tr); 


  partitionList *partitions =  (partitionList*)exa_calloc(1, sizeof(partitionList));
  makeFileNames();

  initializeTree(tr, partitions, adef); 
  printModelAndProgramInfo(tr, partitions, adef, argc, argv);
  exa_main(tr, partitions, adef); 

  return 0;
}


