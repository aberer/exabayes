#include "config.h"
#include "mem_alloc.h"
#include "common.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"
#undef GLOBAL_VARIABLES_DEFINITION

#include "main-common.h"

#include "proposals.h"
#include "output.h"
/* turn on, when in release mode */
/* #define PRODUCTIVE */



/* TODO this is really bad... we really have to clean this up... */
/* more global variables =(  */
char configFileName[1024]; 
int seed; 
int Thorough = 0;
int processID;

void mcmc(tree *tr, analdef *adef); 
boolean setupTree (tree *tr, boolean doInit);
void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile); 


/* TODO: Modify this to our needs */
static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{

  int i, model;
  char modelType[128];


  if(tr->useMedian)
    strcpy(modelType, "GAMMA with Median");
  else
    strcpy(modelType, "GAMMA");

  printVersionInfo();

  if(!adef->compressPatterns)
    PRINT( "\nAlignment has %d columns\n\n",  tr->originalCrunchedLength);
  else
    PRINT( "\nAlignment has %d distinct alignment patterns\n\n",  tr->originalCrunchedLength);



  PRINT( "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");


  switch(adef->mode)
  {	
    case  BIG_RAPID_MODE:	 
      PRINT("\nRAxML rapid hill-climbing mode\n\n");
      break;	
    case  GPU_BENCHMARK:	 
      PRINT("\nRAxML GPU benchmark\n\n");
      break;	
    default:
      assert(0);
  }

  if(adef->perGeneBranchLengths)
    PRINT( "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
  else
    PRINT( "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);	

  PRINT( "All free model parameters will be estimated by RAxML\n");

  if(tr->rateHetModel == GAMMA)
    PRINT("%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
  else    
    PRINT("ML estimate of %d per site rate categories\n\n", tr->categories);

  for(model = 0; model < tr->NumberOfModels; model++)
  {
    PRINT("Partition: %d\n", model);
    PRINT("Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
    PRINT("Name: %s\n", tr->partitionData[model].partitionName);

    switch(tr->partitionData[model].dataType)
    {
      case DNA_DATA:
        PRINT("DataType: DNA\n");	     
        PRINT("Substitution Matrix: GTR\n");
        break;
      case AA_DATA:
        assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
        PRINT("DataType: AA\n");	      
        PRINT("Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);
        PRINT("%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");	     
        break;
      case BINARY_DATA:
        PRINT("DataType: BINARY/MORPHOLOGICAL\n");	      
        PRINT("Substitution Matrix: Uncorrected\n");
        break;
      case SECONDARY_DATA:
        PRINT("DataType: SECONDARY STRUCTURE\n");	     
        PRINT("Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
        break;
      case SECONDARY_DATA_6:
        PRINT("DataType: SECONDARY STRUCTURE 6 STATE\n");	     
        PRINT("Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
        break;
      case SECONDARY_DATA_7:
        PRINT("DataType: SECONDARY STRUCTURE 7 STATE\n");	      
        PRINT("Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
        break;
      case GENERIC_32:
        PRINT("DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
        switch(tr->multiStateModel)
        {
          case ORDERED_MULTI_STATE:
            PRINT("Substitution Matrix: Ordered Likelihood\n");
            break;
          case MK_MULTI_STATE:
            PRINT("Substitution Matrix: MK model\n");
            break;
          case GTR_MULTI_STATE:
            PRINT("Substitution Matrix: GTR\n");
            break;
          default:
            assert(0);
        }
        break;
      case GENERIC_64:
        PRINT("DataType: Codon\n");		  
        break;		
      default:
        assert(0);
    }
    PRINT("\n\n\n");
  }

  PRINT("\n");

  PRINT("%s was called as follows:\n\n", PROGRAM_NAME);
  for(i = 0; i < argc; i++)
    PRINT("%s ", argv[i]);
  PRINT("\n\n\n");

}


int main (int argc, char *argv[])
{ 
    
  tree  *tr = (tree*)rax_malloc(sizeof(tree));

  analdef *adef = (analdef*)rax_malloc(sizeof(analdef));

  double **empiricalFrequencies;

  ignoreExceptionsDenormFloat();

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
  if(workerTrap(tr))
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

  
  size_t 
    i,
    model;

  unsigned char *y;

  FILE 
    *byteFile = myfopen(byteFileName, "rb");	 

  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
  myBinFread(&(tr->NumberOfModels),         sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);

  if(adef->perGeneBranchLengths)
    tr->numBranches = tr->NumberOfModels;
  else
    tr->numBranches = 1;

  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed RAxML-specific */

  tr->aliaswgt                   = (int *)rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
  myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       

  tr->rateCategory    = (int *)    rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));	  

  tr->patrat          = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->patratStored    = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
  tr->lhs             = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 

  tr->executeModel   = (boolean *)rax_malloc(sizeof(boolean) * (size_t)tr->NumberOfModels);

  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    tr->executeModel[i] = TRUE;

  empiricalFrequencies = (double **)rax_malloc(sizeof(double *) * (size_t)tr->NumberOfModels);

  y = (unsigned char *)rax_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));

  tr->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength]; 

  setupTree(tr, FALSE);

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

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      int 
        len;

      pInfo 
        *p = &(tr->partitionData[model]);	   

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

      empiricalFrequencies[model] = (double *)rax_malloc(sizeof(double) * (size_t)tr->partitionData[model].states);
      myBinFread(empiricalFrequencies[model], sizeof(double), tr->partitionData[model].states, byteFile);	   
    }

  myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

  fclose(byteFile);


  /* Now we are able to compute the memory requirements and decide on using recom or not */
  if(tr->useRecom)
    {
      size_t requiredLength;
      size_t rateHet = discreteRateCategories(tr->rateHetModel);
      int model;
      float approxTotalMegabytesRequired;
      requiredLength = 0;
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t width  = (size_t)tr->partitionData[model].width;
	  size_t states = (size_t)tr->partitionData[model].states;
	  requiredLength += virtual_width(width) * rateHet * states * sizeof(double);
	}
      requiredLength *= tr->mxtips - 2;
      approxTotalMegabytesRequired = MEM_APROX_OVERHEAD * (requiredLength / 1E6);  
      PRINT("Required memory for inner nodes in MB: %f\n", (float)requiredLength / 1E6);
      PRINT("Estimated total required memory in MB: %f\n", approxTotalMegabytesRequired);
      PRINT("User      maximum use of memory in MB: %f\n", tr->maxMegabytesMemory);
      if (approxTotalMegabytesRequired < tr->maxMegabytesMemory)
	{
	  tr->useRecom = FALSE; 
	  PRINT("Deactivated recomputation of inner vectors\n");
	}
      else
	{
	  tr->vectorRecomFraction = tr->maxMegabytesMemory  / approxTotalMegabytesRequired;
	  PRINT("Set recomputation of inner vectors to fraction %f\n", tr->vectorRecomFraction);
	  assert(tr->vectorRecomFraction > MIN_RECOM_FRACTION);
	  assert(tr->vectorRecomFraction < MAX_RECOM_FRACTION);
	}
    }


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* 
     this main function is the master thread, so if we want to run RAxML with n threads,
     we use startPthreads to start the n-1 worker threads */
  
#ifdef _USE_PTHREADS
  startPthreads(tr);
#endif

  /* via masterBarrier() we invoke parallel regions in which all Pthreads work on computing something, mostly likelihood 
     computations. Have a look at execFunction() in axml.c where we siwtch of the different types of parallel regions.

     Although not necessary, below we copy the info stored on tr->partitionData to corresponding copies in each thread.
     While this is shared memory and we don't really need to copy stuff, it was implemented like this to allow for an easier 
     transition to a distributed memory implementation (MPI).
  */
  
  /* mpi version now also uses the generic barrier */
  masterBarrier(THREAD_INIT_PARTITION, tr); 
#else  /* SEQUENTIAL */
  /* 
     allocate the required data structures for storing likelihood vectors etc 
  */

  initializePartitions(tr, tr, 0, 0);
#endif

  /* print out some info on partitions, models, data types etc, not very interesting */

  printModelAndProgramInfo(tr, adef, argc, argv);

  /* Tells us if the SEV-based memory saving option has been activated in the command line or not.
     PRINT() allows to simultaneously print to terminal and to the RAxML_info file, thereby making 
     the terminal output and the info in the RAxML_info file consistent */

  PRINT("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             

  /* Tells us if the vector recomputation memory saving option has been activated in the command line or not.
   */
  /* recom */
  
  PRINT("Memory Saving via Additional Vector Recomputations: %s\n", (tr->useRecom == TRUE)?"ENABLED":"DISABLED");
  if(tr->useRecom)
    PRINT("Using a fraction %f of the total inner vectors that would normally be required\n", tr->vectorRecomFraction);
    
  /* E recom */

  /* 
     initialize model parameters like empirical base frequencies, the rates in the Q matrix, the alpha shape parameters,
     the per-site substitution rates to default starting values */

  initModel(tr, empiricalFrequencies); 
  PRINT("Model initialized\n");

  /* 
     this will re-start RAxML exactly where it has left off from a checkpoint file,
     while checkpointing is important and has to be implemented for the library we should not worry about this right now 
  */
  

  /* not important, only used to keep track of total accumulated exec time 
     when checkpointing and restarts were used */

  accumulatedTime = 0.0;

  /* get the starting tree: here we just parse the tree passed via the command line 
     and do an initial likelihood computation traversal 
     which we maybe should skeip, TODO */
  
  
  tr->startingTree = 1 ; 
  getStartingTree(tr);

  /* 
     here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
     This should basically be the first call to the library that actually computes something :-)
  */


  /* please do not remove this code from here ! */

  evaluateGeneric(tr, tr->start, TRUE);
  PRINT("Starting tree evaluated\n");

  /* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */

  treeEvaluate(tr, 32);
  PRINT("tree evaluated: %f\n", tr->likelihood);
    
  /* now start the ML search algorithm */

  /* allocate parsimony data structures for parsimony-biased SPRs */	
  allocateParsimonyDataStructures(tr);

  mcmc(tr, adef);

  freeParsimonyDataStructures(tr);

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* workers escape from their while loop (could be joined in pthread case )  */
  masterBarrier(THREAD_EXIT_GRACEFULLY,tr);
#endif

  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  return 0;
}

