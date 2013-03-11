#include "config.h"
#include "mem_alloc.h"
#include "common.h"

#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"
#undef GLOBAL_VARIABLES_DEFINITION

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "main-common.h"

#include "proposalStructs.h"
#include "proposals.h"
#include "output.h"
/* turn on, when in release mode */
/* #define PRODUCTIVE */


#include "adapterCode.h"

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
  tree  *tr = (tree*)rax_malloc(sizeof(tree));
  /* partitionList *partitions = (partitionList*)rax_malloc(sizeof(partitionList)); */
  partitionList *partitions = &(gAInfo.partitions); 

  partitions->partitionData = (pInfo**)rax_malloc(NUM_BRANCHES*sizeof(pInfo*));

  analdef *adef = (analdef*)rax_malloc(sizeof(analdef));

  double **empiricalFrequencies;

  /* 
     tell the CPU to ignore exceptions generated by denormalized floating point values.
     If this is not done, depending on the input data, the likelihood functions can exhibit 
     substantial run-time differences for vectors of equal length.
  */

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
  _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif 


  /* get the start time */

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

    tr->aliaswgt                   = (int *)rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
    myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       

    tr->rateCategory    = (int *)    rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));	  

    tr->patrat          = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
    tr->patratStored    = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
    tr->lhs             = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 

    empiricalFrequencies = (double **)rax_malloc(sizeof(double *) * partitions->numberOfPartitions);

    y = (unsigned char *)rax_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));
    tr->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

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
	tr->nameList[i] = (char*)rax_malloc(sizeof(char) * (size_t)len);
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
	p->partitionName = (char*)rax_malloc(sizeof(char) * (size_t)len);
	myBinFread(p->partitionName, sizeof(char), len, byteFile);

	empiricalFrequencies[model] = (double *)rax_malloc(sizeof(double) * (size_t)partitions->partitionData[model]->states);
	myBinFread(empiricalFrequencies[model], sizeof(double), partitions->partitionData[model]->states, byteFile);
      }

    myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

    fclose(byteFile);
  }


  /* gAInfo.partitions = partitions; */


  /* Now we are able to compute the memory requirements and decide on using recom or not */
  if(tr->useRecom)
    {
      size_t requiredLength;
      size_t rateHet = discreteRateCategories(tr->rateHetModel);
      int model;
      float approxTotalMegabytesRequired;
      requiredLength = 0;
      for(model = 0; model < partitions->numberOfPartitions; model++)
	{
	  size_t width  = (size_t)partitions->partitionData[model]->width;
	  size_t states = (size_t)partitions->partitionData[model]->states;
	  requiredLength += virtual_width(width) * rateHet * states * sizeof(double);
	}
      requiredLength *= tr->mxtips - 2;
      approxTotalMegabytesRequired = MEM_APROX_OVERHEAD * (requiredLength / 1E6);  
      printBothOpen("Required memory for inner nodes in MB: %f\n", (float)requiredLength / 1E6);
      printBothOpen("Estimated total required memory in MB: %f\n", approxTotalMegabytesRequired);
      printBothOpen("User      maximum use of memory in MB: %f\n", tr->maxMegabytesMemory);
      if (approxTotalMegabytesRequired < tr->maxMegabytesMemory)
	{
	  tr->useRecom = FALSE; 
	  printBothOpen("Deactivated recomputation of inner vectors\n");
	}
      else
	{
	  tr->vectorRecomFraction = tr->maxMegabytesMemory  / approxTotalMegabytesRequired;
	  printBothOpen("Set recomputation of inner vectors to fraction %f\n", tr->vectorRecomFraction);
	  assert(tr->vectorRecomFraction > MIN_RECOM_FRACTION);
	  assert(tr->vectorRecomFraction < MAX_RECOM_FRACTION);
	}
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

  /* Tells us if the SEV-based memory saving option has been activated in the command line or not.
     printBothOpen() allows to simultaneously print to terminal and to the RAxML_info file, thereby making 
     the terminal output and the info in the RAxML_info file consistent */

  printBothOpen("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             

  /* Tells us if the vector recomputation memory saving option has been activated in the command line or not.
   */
  /* recom */
  {
#ifdef _DEBUG_RECOMPUTATION
    allocTraversalCounter(tr);
    tr->stlenTime = 0.0;
    if(tr->useRecom)
      {
	tr->rvec->pinTime = 0.0;
	tr->rvec->recomStraTime = 0.0;
      }
#endif
    printBothOpen("Memory Saving via Additional Vector Recomputations: %s\n", (tr->useRecom == TRUE)?"ENABLED":"DISABLED");
    if(tr->useRecom)
      printBothOpen("Using a fraction %f of the total inner vectors that would normally be required\n", tr->vectorRecomFraction);
  }
  /* E recom */

  /* 
     initialize model parameters like empirical base frequencies, the rates in the Q matrix, the alpha shape parameters,
     the per-site substitution rates to default starting values */

  initModel(tr, empiricalFrequencies, partitions);
  printBothOpen("Model initialized\n");
  




  /* 
     this will re-start RAxML exactly where it has left off from a checkpoint file,
     while checkpointing is important and has to be implemented for the library we should not worry about this right now 
  */

  accumulatedTime = 0.0;

  /* get the starting tree: here we just parse the tree passed via the command line 
     and do an initial likelihood computation traversal 
     which we maybe should skeip, TODO */

  printBothOpen("Starting tree available\n");
  /* 
     here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
     This should basically be the first call to the library that actually computes something :-)
  */


  /* please do not remove this code from here ! */

  /* evaluateGeneric(tr, partitions, tr->start, TRUE); */
  /* printBothOpen("Starting tree evaluated\n"); */


  /**** test code for testing per-site log likelihood calculations as implemented in evaluatePartialGenericSpecial.c for Kassian's work*/

 
  /* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */

  /* treeEvaluate(tr, partitions, 32); */
  /* printBothOpen("tree evaluated: %f\n", tr->likelihood); */
    
  /* now start the ML search algorithm */
  printf("starting MCMC\n"); 
  
  mcmc(tr, adef);

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


