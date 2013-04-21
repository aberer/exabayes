#ifdef HAVE_AVX
#define __AVX
#endif

#include "axml.h" 
#include "bayes.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "main-common.h"
#include "proposals.h"
#include "output.h"
#include "adapters.h"

#if HAVE_PLL != 0

#include "globalVariables.h"

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
#if HAVE_PLL != 0
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




#else 


extern int processID; 
extern char infoFileName[1024];
extern char *protModels[NUM_PROT_MODELS]; 
extern const char *secondaryModelList[21]; 
extern partitionLengths pLengths[MAX_MODEL]; 
extern partitionLengths pLength;
extern char resultFileName[1024]; 
extern double masterTime;
extern char logFileName[1024]; 
extern int processes;
extern double accumulatedTime;
extern char byteFileName[1024]; 

void exa_main(tree *tr, analdef *adef); 



void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[]);
void initializeTree(tree *tr, analdef *adef);


int main(int argc, char *argv[])
{ 
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  printf("\nThis is %s FINE-GRAIN MPI Process Number: %d\n", PROGRAM_NAME, processID);   
  MPI_Barrier(MPI_COMM_WORLD);

  tree  *tr = (tree*)exa_malloc(sizeof(tree));  
  analdef *adef = (analdef*)exa_malloc(sizeof(analdef));   
  
  ignoreExceptionsDenormFloat(); 

  /* get the start time */
   
  masterTime = gettime();         
    
  /* initialize the analysis parameters in struct adef to default values */
    
  initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */
  
  get_args(argc, argv, adef, tr); 
  
  /* generate the ExaML output file names and store them in strings */
    
  makeFileNames();
  
  initializeTree(tr, adef); 

  if(isOutputProcess()) 
    {
      printModelAndProgramInfo(tr, adef, argc, argv);
      PRINT("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             
    }  
                         
  /* 
     this will re-start ExaML exactly where it has left off from a checkpoint file,
     while checkpointing is important and has to be implemented for the library we should not worry about this right now 
  */

  /* not important, only used to keep track of total accumulated exec time 
     when checkpointing and restarts were used */
	
  if(isOutputProcess())
    accumulatedTime = 0.0;
 
  exa_main(tr,adef);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}


#endif




