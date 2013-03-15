#include "config.h"
#include "mem_alloc.h"
#include "common.h"

#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"
#undef GLOBAL_VARIABLES_DEFINITION

#include "proposalStructs.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "main-common.h"


#include "proposals.h"
#include "output.h"
/* turn on, when in release mode */
/* #define PRODUCTIVE */


#include "adapterCode.h"

void exa_main(tree *tr, analdef *adef); 

boolean setupTree (tree *tr, boolean doInit, partitionList *partitions);

/* void initializePartitions(tree *tr, tree *localTree, partitionList *pr, partitionList *localPr, int tid, int n);  */

/* TODO this is really bad... we really have to clean this up... */
/* more global variables =(  */
char configFileName[1024]; 
int seed; 
int Thorough = 0;
int processID;

void mcmc(tree *tr, analdef *adef); 
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

static void finalizeInfoFile(tree *tr, analdef *adef)
{
  double t;

  t = gettime() - masterTime;
  accumulatedTime = accumulatedTime + t;

  switch(adef->mode)
  {
    case  BIG_RAPID_MODE:
      printBothOpen("\n\nOverall Time for 1 Inference %f\n", t);
      printBothOpen("\nOverall accumulated Time (in case of restarts): %f\n\n", accumulatedTime);
      printBothOpen("Likelihood   : %f\n", tr->likelihood);
      printBothOpen("\n\n");
      printBothOpen("Final tree written to:                 %s\n", resultFileName);
      printBothOpen("Execution Log File written to:         %s\n", logFileName);
      printBothOpen("Execution information file written to: %s\n",infoFileName);
      break;
    default:
      assert(0);
  }
}



static void printModelAndProgramInfo(tree *tr, partitionList *pr, analdef *adef, int argc, char *argv[])
{

  int i, model;
  FILE *infoFile = myfopen(infoFileName, "ab");
  char modelType[128];


  if(tr->useMedian)
    strcpy(modelType, "GAMMA with Median");
  else
    strcpy(modelType, "GAMMA");

  printBoth(infoFile, "\n\nThis is %s version %s released by Alexandros Stamatakis in %s.\n\n",  programName, programVersion, programDate);



  if(!adef->compressPatterns)
    printBoth(infoFile, "\nAlignment has %d columns\n\n",  tr->originalCrunchedLength);
  else
    printBoth(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->originalCrunchedLength);



  printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");


  switch(adef->mode)
  {	
    case  BIG_RAPID_MODE:	 
      printBoth(infoFile, "\nRAxML rapid hill-climbing mode\n\n");
      break;	
    case  GPU_BENCHMARK:	 
      printBoth(infoFile, "\nRAxML GPU benchmark\n\n");
      break;	
    default:
      assert(0);
  }

  if(adef->perGeneBranchLengths)
    printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", pr->numberOfPartitions);
  else
    printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", pr->numberOfPartitions);

  printBoth(infoFile, "All free model parameters will be estimated by RAxML\n");

  if(tr->rateHetModel == GAMMA)
    printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
  else    
    printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", tr->categories);

  for(model = 0; model < pr->numberOfPartitions; model++)
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

  printBoth(infoFile, "RAxML was called as follows:\n\n");
  for(i = 0; i < argc; i++)
    printBoth(infoFile,"%s ", argv[i]);
  printBoth(infoFile,"\n\n\n");

  fclose(infoFile);
}



int main (int argc, char *argv[])
{ 
  tree  *tr = (tree*)exa_malloc(sizeof(tree));
  /* partitionList *partitions = (partitionList*)exa_malloc(sizeof(partitionList)); */
  partitionList *partitions = gAInfo.partitions; 

  partitions->partitionData = (pInfo**)exa_malloc(NUM_BRANCHES*sizeof(pInfo*));

  analdef *adef = (analdef*)exa_malloc(sizeof(analdef));

  double **empiricalFrequencies;
  ignoreExceptionsDenormFloat(); 

  masterTime = gettime();         

  /* 
     initialize the workers for mpi or pthreads
  */
#ifdef _FINE_GRAIN_MPI
  /* 
     once mpi workers are signalled to finish, it is impontant that
     they immediately terminate! (to avoid undefined behavior)
  */
#ifdef MEASURE_TIME_PARALLEL
  masterTimePerPhase = gettime();
#endif
  initMPI(argc, argv);
  if(workerTrap(tr, partitions))
    return 0; 
#endif
#ifdef _USE_PTHREADS
  tr->threadID = 0;
#ifndef _PORTABLE_PTHREADS
  /* not very portable thread to core pinning if PORTABLE_PTHREADS is not defined
     by defualt the cod ebelow is deactivated */
  pinToCore(0);
#endif
#endif


  /* initialize the analysis parameters in struct adef to default values */

  initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */

  get_args(argc, argv, adef, tr); 

  /* generate the RAxML output file names and store them in strings */

  makeFileNames();

  {
    size_t 
      i,
      model;

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

    empiricalFrequencies = (double **)exa_malloc(sizeof(double *) * partitions->numberOfPartitions);

    y = (unsigned char *)exa_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));
    tr->yVector = (unsigned char **)exa_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength]; 

    setupTree(tr, FALSE, partitions);

    for(i = 0; i < partitions->numberOfPartitions; i++)
      partitions->partitionData[i]->executeModel = TRUE;

    /* data structures for convergence criterion need to be initialized after! setupTree */

    if(tr->searchConvergenceCriterion)
      {                     
	tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
	tr->h = initHashTable(tr->mxtips * 4);        
      }

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      {
	int len;
	myBinFread(&len, sizeof(int), 1, byteFile);
	tr->nameList[i] = (char*)exa_malloc(sizeof(char) * (size_t)len);
	myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
	/*printf("%s \n", tr->nameList[i]);*/
      }  

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      addword(tr->nameList[i], tr->nameHash, i);

    for(model = 0; model < partitions->numberOfPartitions; model++)
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

	empiricalFrequencies[model] = (double *)exa_malloc(sizeof(double) * (size_t)partitions->partitionData[model]->states);
	myBinFread(empiricalFrequencies[model], sizeof(double), partitions->partitionData[model]->states, byteFile);
      }

    myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

    fclose(byteFile);
  }

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* 
     this main function is the master thread, so if we want to run RAxML with n threads,
     we use startPthreads to start the n-1 worker threads */
  
#ifdef _USE_PTHREADS
  startPthreads(tr, partitions);
#endif

  /* via masterBarrier() we invoke parallel regions in which all Pthreads work on computing something, mostly likelihood 
     computations. Have a look at execFunction() in axml.c where we siwtch of the different types of parallel regions.

     Although not necessary, below we copy the info stored on tr->partitionData to corresponding copies in each thread.
     While this is shared memory and we don't really need to copy stuff, it was implemented like this to allow for an easier 
     transition to a distributed memory implementation (MPI).
  */
  
  /* mpi version now also uses the generic barrier */
  masterBarrier(THREAD_INIT_PARTITION, tr, partitions);
#else  /* SEQUENTIAL */
  /* 
     allocate the required data structures for storing likelihood vectors etc 
  */

  initializePartitions(tr, tr, partitions, 0, 0);
#endif

  /* print out some info on partitions, models, data types etc, not very interesting */

  printModelAndProgramInfo(tr, partitions, adef, argc, argv);

  initModel(tr, empiricalFrequencies, partitions);

  exa_main(tr, adef); 
  


#ifdef _DEBUG_RECOMPUTATION
  {
    double t = gettime() - masterTime;
    printBothOpen("Traversal freq after search \n");
    printTraversalInfo(tr);
    if(tr->useRecom)
      printBothOpen("Recom stlen %f, cost %f, pin %f, t %f\n", 
		    tr->stlenTime,tr->rvec->recomStraTime, tr->rvec->pinTime, t);
    else
      printBothOpen("No Recom stlen %f, t %f\n", tr->stlenTime, t);
  }
#endif 

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* workers escape from their while loop (could be joined in pthread case )  */
  masterBarrier(THREAD_EXIT_GRACEFULLY,tr, partitions);
#endif

  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  return 0;
}


