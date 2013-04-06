

/* use, if in productive mode */
/* #define PRODUCTIVE */


#include <mpi.h>

#include "axml.h"
#include "globalVariables.h"
#include "main-common.h"

#include "bayes.h"
#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "bayes.h"
#include "output.h"

#include "adapters.h"


void exa_main(tree *tr, analdef *adef); 
void initAdef(analdef *adef);
void makeFileNames(void); 
void initializeTree(tree *tr, analdef *adef); 
void finalizeInfoFile(tree *tr, analdef *adef); 


void makeRandomTree(tree *tr); 



/***************** UTILITY FUNCTIONS **************************/

void storeExecuteMaskInTraversalDescriptor(tree *tr)
{
   int model;
      
   for(model = 0; model < getNumberOfPartitions(tr); model++)
     tr->td[0].executeModel[model] = tr->executeModel[model];
}

void storeValuesInTraversalDescriptor(tree *tr, double *value)
{
   int model;
      
   for(model = 0; model < getNumberOfPartitions(tr); model++)
     tr->td[0].parameterValues[model] = value[model];
}



void printBothOpen(const char* format, ... )
{
  if(processID == 0)
    {
      FILE *f = myfopen(infoFileName, "ab");
      
      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
      
      fclose(f);
    }
}

static void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{
  size_t
    bytes_read;
  
  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{
  int numberOfPartitions = getNumberOfPartitions(tr); 

  if(processID == 0)
    {
      int i, model;
      FILE *infoFile = myfopen(infoFileName, "ab");
      char modelType[128];

      
      if(tr->useMedian)
	strcpy(modelType, "GAMMA with Median");
      else
	strcpy(modelType, "GAMMA");   

      PRINT("\n\nThis is %s version %s.\n\n", PROGRAM_NAME, PACKAGE_VERSION); 
      
      PRINT( "\nAlignment has %z distinct alignment patterns\n\n",  tr->originalCrunchedLength);
      
      PRINT( "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");


      if(adef->perGeneBranchLengths)
	PRINT( "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", numberOfPartitions);
      else
	PRINT( "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", numberOfPartitions);	


      /*
	if(adef->mode != CLASSIFY_ML)
	PRINT( "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n",
	modelType, adef->likelihoodEpsilon);
      */
    
      
      for(model = 0; model < numberOfPartitions; model++)
	{
	  PRINT( "Partition: %d\n", model);
	  PRINT( "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
	  PRINT( "Name: %s\n", tr->partitionData[model].partitionName);
	  
	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      PRINT( "DataType: DNA\n");	     
	      PRINT( "Substitution Matrix: GTR\n");
	      break;
	    case AA_DATA:
	      assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
	      PRINT( "DataType: AA\n");	      
	      PRINT( "Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);
	      PRINT( "%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");	     
	      break;
	    case BINARY_DATA:
	      PRINT( "DataType: BINARY/MORPHOLOGICAL\n");	      
	      PRINT( "Substitution Matrix: Uncorrected\n");
	      break;
	    case SECONDARY_DATA:
	      PRINT( "DataType: SECONDARY STRUCTURE\n");	     
	      PRINT( "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	      break;
	    case SECONDARY_DATA_6:
	      PRINT( "DataType: SECONDARY STRUCTURE 6 STATE\n");	     
	      PRINT( "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	      break;
	    case SECONDARY_DATA_7:
	      PRINT( "DataType: SECONDARY STRUCTURE 7 STATE\n");	      
	      PRINT( "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	      break;
	    case GENERIC_32:
	      PRINT( "DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
	      switch(tr->multiStateModel)
		{
		case ORDERED_MULTI_STATE:
		  PRINT( "Substitution Matrix: Ordered Likelihood\n");
		  break;
		case MK_MULTI_STATE:
		  PRINT( "Substitution Matrix: MK model\n");
		  break;
		case GTR_MULTI_STATE:
		  PRINT( "Substitution Matrix: GTR\n");
		  break;
		default:
		  assert(0);
		}
	      break;
	    case GENERIC_64:
	      PRINT( "DataType: Codon\n");		  
	      break;		
	    default:
	      assert(0);
	    }
	  PRINT( "\n\n\n");
	}
      
      PRINT( "\n");

      PRINT( "%s was called as follows:\n\n", PROGRAM_NAME);
      for(i = 0; i < argc; i++)
	PRINT("%s ", argv[i]);
      PRINT("\n\n\n");

      fclose(infoFile);
    }
}

void *malloc_aligned(size_t size) 
{
  void 
    *ptr = (void *)NULL;
 
  int 
    res;
  

#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
  */

  ptr = exa_malloc(size);
  
  if(ptr == (void*)NULL) 
   assert(0);
  
#ifdef __AVX
  assert(0);
#endif

#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0) 
    assert(0);
#endif 

  return ptr;
}


boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}



char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}

partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;
  pLength.nonGTR = FALSE;

  return (&pLengths[dataType]); 
}










size_t discreteRateCategories(int rateHetModel)
{
  size_t 
    result;

  switch(rateHetModel)
    {
    case CAT:
      result = 1;
      break;
    case GAMMA:
      result = 4;
      break;
    default:
      assert(0);
    }

  return result;
}



double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

/* int gettimeSrand(void) */
/* { */
/* #ifdef WIN32 */
/*   time_t tp; */
/*   struct tm localtm; */
/*   tp = time(NULL); */
/*   localtm = *localtime(&tp); */
/*   return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec; */
/* #else */
/*   struct timeval ttime; */
/*   gettimeofday(&ttime , NULL); */
/*   return ttime.tv_sec + ttime.tv_usec; */
/* #endif */
/* } */

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s %s wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, PROGRAM_NAME, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }

}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}









void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }

  assert(p->x);
}





void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}


/***********************reading and initializing input ******************/







boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}










static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}

static boolean setupTree (tree *tr)
{
  nodeptr  p0, p, q;
  int
    i,
    j,   
    tips,
    inter; 

  int numPartitions = getNumberOfPartitions(tr) ;

  
  tr->bigCutoff = FALSE;
  
  tr->maxCategories = MAX(4, tr->categories);
  
  tr->partitionContributions = (double *)exa_malloc(sizeof(double) * numPartitions);
  
  for(i = 0; i < getNumberOfPartitions(tr); i++)
    tr->partitionContributions[i] = -1.0;
  
  tr->perPartitionLH = (double *)exa_malloc(sizeof(double) * numPartitions);
  
  
  for(i = 0; i < numPartitions; i++)    
    tr->perPartitionLH[i] = 0.0;	    
  
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  tr->fracchanges  = (double *)exa_malloc(getNumberOfPartitions(tr) * sizeof(double));
  

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)exa_calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)exa_calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)exa_calloc(tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  
            
  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)exa_malloc(sizeof(traversalInfo) * tr->mxtips);
  tr->td[0].executeModel = (boolean *)exa_malloc(sizeof(boolean) * getNumberOfPartitions(tr));
  tr->td[0].parameterValues = (double *)exa_malloc(sizeof(double) * getNumberOfPartitions(tr));
  
  for(i = 0; i < getNumberOfPartitions(tr); i++)
    tr->fracchanges[i] = -1.0;
  tr->fracchange = -1.0;
  
  tr->constraintVector = (int *)exa_malloc((2 * tr->mxtips) * sizeof(int));
  
  tr->nameList = (char **)exa_malloc(sizeof(char *) * (tips + 1));
   

  if (!(p0 = (nodeptr) exa_malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }
  
  tr->nodeBaseAddress = p0;


  if (!(tr->nodep = (nodeptr *) exa_malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;
      p->xBips  =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;     
      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    {
	      p->xBips = 1;
	      p->x = 1;
	    }
	  else
	    {
	      p->xBips = 0;
	      p->x =  0;
	    }
	  p->number = i;
	  p->next   = q;	  
	  p->back   = (node *) NULL;
	  p->hash   = 0;       
	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;
  
  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  tr->partitionData = (pInfo*)exa_malloc(sizeof(pInfo) * getNumberOfPartitions(tr));

  return TRUE;
}



/*********************************** *********************************************************/




/********************PRINTING various INFO **************************************/


void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0)
    {
      FILE *logFile;
      char temporaryFileName[1024] = "";
      
      strcpy(temporaryFileName, resultFileName);
      
      switch(adef->mode)
	{    
	case TREE_EVALUATION:
	  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, SUMMARIZE_LH, FALSE, FALSE);
	  
	  logFile = myfopen(temporaryFileName, "wb");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	  
	  if(adef->perGeneBranchLengths)
	    printTreePerGene(tr, adef, temporaryFileName, "wb");
	  break;
	case BIG_RAPID_MODE:     
	  if(finalPrint)
	    {
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		case GAMMA_I:
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint,
			      SUMMARIZE_LH, FALSE, FALSE);
		  
		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "wb");
		  break;
		case CAT:
		  /*Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
		    NO_BRANCoHES, FALSE, FALSE);*/
		  
		  
		  
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE,
			      TRUE, SUMMARIZE_LH, FALSE, FALSE);
		  
		  
		  
		  
		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  break;
		default:
		  assert(0);
		}
	    }
	  else
	    {
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint,
			  NO_BRANCHES, FALSE, FALSE);
	      logFile = myfopen(temporaryFileName, "wb");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);
	    }    
	  break;
	default:
	  printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
	  exit(-1);
	  break;
	}
    }
}








void printLog(tree *tr)
{
  if(processID == 0)
    {
      FILE *logFile;
      double t;
      
      t = gettime() - masterTime;
      
      logFile = myfopen(logFileName, "ab");
      
      /* printf("%f %1.40f\n", t, tr->likelihood); */

      fprintf(logFile, "%f %f\n", t, tr->likelihood);
      
      fclose(logFile);
    }
	     
}









void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
    {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
    }
}

 


void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(processID == 0)
    {
      double t;

      t = gettime() - masterTime;
      accumulatedTime = accumulatedTime + t;

      switch(adef->mode)
	{	
	case  BIG_RAPID_MODE:	 
	  PRINT("\n\nOverall Time for 1 Inference %f\n", t);
	  PRINT("\nOverall accumulated Time (in case of restarts): %f\n\n", accumulatedTime);
	  PRINT("Likelihood   : %f\n", tr->likelihood);
	  PRINT("\n\n");	  	  
	  PRINT("Final tree written to:                 %s\n", resultFileName);
	  PRINT("Execution Log File written to:         %s\n", logFileName);
	  PRINT("Execution information file written to: %s\n",infoFileName);	
	  break;
	default:
	  assert(0);
	}

	 
    }

}


/************************************************************************************/







boolean isThisMyPartition(tree *tr, int tid, int model)
{ 
  assert(tr->manyPartitions);

  if(tr->partitionAssignment[model] == tid)
    return TRUE;
  else
    return FALSE;
}

static void computeFractionMany(tree *tr, int tid)
{
  int
    sites = 0;

  int   
    model;

  assert(tr->manyPartitions);

  for(model = 0; model < getNumberOfPartitions(tr); model++)
    {
      if(isThisMyPartition(tr, tid, model))
	{	 
	  tr->partitionData[model].width = tr->partitionData[model].upper - tr->partitionData[model].lower;
	  sites += tr->partitionData[model].width;
	}
      else       	  
	tr->partitionData[model].width = 0;       
    }

  
}



static void computeFraction(tree *tr, int tid, int n)
{
  int
    model;

  size_t 
    i;

  for(model = 0; model < getNumberOfPartitions(tr); model++)
    {
      size_t 
	width = 0;

      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	if(i % n == (size_t)tid)
	  width++;

      tr->partitionData[model].width = width;
    }
}








static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

/*static char bits_in_16bits [0x1u << 16];*/

static void compute_bits_in_16bits(char *bits_in_16bits)
{
    unsigned int i;    
    
    /* size is 65536 */

    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);       

    return ;
}

unsigned int precomputed16_bitcount (unsigned int n, char *bits_in_16bits)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}



 

static int partCompare(const void *p1, const void *p2)
{
 partitionType 
   *rc1 = (partitionType *)p1,
   *rc2 = (partitionType *)p2;

 int 
   i = rc1->partitionLength,
   j = rc2->partitionLength;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


/* 
	 tr->manyPartitions is set to TRUE if the user has indicated via -Q that there are substantially more partitions 
	 than threads/cores available. In that case we do not distribute sites from each partition in a cyclic fashion to the cores 
	 , but distribute entire partitions to cores. 
	 Achieving a good balance of alignment sites to cores boils down to the mult-processor scheduling problem known from theoretical comp. sci.
	 which si NP-complete.
	 We have implemented very simple "standard" heuristics for solving the multiprocessor scheduling problem that turn out to work very well
	 and are cheap to compute 
*/

static void multiprocessorScheduling(tree *tr, int tid)
{
  int 
    s,
    model,
    modelStates[2] = {4, 20},
    numberOfPartitions[2] = {0 , 0},
    arrayLength = sizeof(modelStates) / sizeof(int);
  
    /* check that we have not addedd any new models for data types with a different number of states
       and forgot to update modelStates */
    
    tr->partitionAssignment = (int *)exa_malloc(getNumberOfPartitions(tr) * sizeof(int));
    
  for(model = 0; model < getNumberOfPartitions(tr); model++)
    {        
      boolean 
	exists = FALSE;

      for(s = 0; s < arrayLength; s++)
	{
	  exists = exists || (tr->partitionData[model].states == modelStates[s]);
	  if(tr->partitionData[model].states == modelStates[s])
	    numberOfPartitions[s] += 1;
	}

      assert(exists);
    }

  if(tid == 0)
    PRINT("\nMulti-processor partition data distribution enabled (-Q option)\n");

  for(s = 0; s < arrayLength; s++)
    {
      if(numberOfPartitions[s] > 0)
	{
	  size_t   
	    checkSum = 0,
	    sum = 0;
	  
	  int    
	    i,
	    k,
	    n = processes,
	    p = numberOfPartitions[s],    
	    *assignments = (int *)exa_calloc(n, sizeof(int));  
	  
	  partitionType 
	    *pt = (partitionType *)exa_malloc(sizeof(partitionType) * p);
	  
	  
	  for(i = 0, k = 0; i < getNumberOfPartitions(tr); i++)
	    {
	      if(tr->partitionData[i].states == modelStates[s])
		{
		  pt[k].partitionNumber = i;
		  pt[k].partitionLength = tr->partitionData[i].upper - tr->partitionData[i].lower;
		  sum += (size_t)pt[k].partitionLength;
		  k++;
		}
	    }
	  
	  assert(k == p);
	  
	  qsort(pt, p, sizeof(partitionType), partCompare);    
	  
	  for(i = 0; i < p; i++)
	    {
	      int 
		k, 
		min = INT_MAX,
		minIndex = -1;
	      
	      for(k = 0; k < n; k++)	
		if(assignments[k] < min)
		  {
		    min = assignments[k];
		    minIndex = k;
		  }
	      
	      assert(minIndex >= 0);
	      
	      assignments[minIndex] +=  pt[i].partitionLength;
	      assert(pt[i].partitionNumber >= 0 && pt[i].partitionNumber < getNumberOfPartitions(tr));
	      tr->partitionAssignment[pt[i].partitionNumber] = minIndex;
	    }
	  
	  if(tid == 0)
	    {
	      for(i = 0; i < n; i++)	       
		PRINT("Process %d has %d sites for %d state model \n", i, assignments[i], modelStates[s]);		  		
	      
	      PRINT("\n");
	    }

	  for(i = 0; i < n; i++)
	    checkSum += (size_t)assignments[i];
	  
	  assert(sum == checkSum);
	  
	  exa_free(assignments);
	  exa_free(pt);
	}
    } 
}



static void initializePartitions(tree *tr, FILE *byteFile)
{ 
  size_t
    i,
    j,
    width,
    model,
    countOffset,
    myLength = 0;

  int
    maxCategories;

  unsigned char 
    *y;

  compute_bits_in_16bits(tr->bits_in_16bits);
  
  for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
    tr->partitionData[model].width        = 0;

  if(tr->manyPartitions)
    {
      multiprocessorScheduling(tr, processID);  
      computeFractionMany(tr, processID);
    }
  else
    computeFraction(tr, processID, processes);
  	   
  maxCategories = tr->maxCategories;

  for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
    {                       
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[model])); 

      width = tr->partitionData[model].width;

      /* tr->partitionData[model].wr = (double *)exa_malloc(sizeof(double) * width); */
      /* tr->partitionData[model].wr2 = (double *)exa_malloc(sizeof(double) * width);      */

     	
      /* 
	 globalScaler needs to be 2 * tr->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> exa_calloc
       */

      tr->partitionData[model].globalScaler    = (unsigned int *)exa_calloc(2 * tr->mxtips, sizeof(unsigned int));  	         

      tr->partitionData[model].left              = (double *)exa_malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].right             = (double *)exa_malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].EIGN              = (double*)exa_malloc(pl->eignLength * sizeof(double));
      tr->partitionData[model].EV                = (double*)exa_malloc_aligned(pl->evLength * sizeof(double));
      tr->partitionData[model].EI                = (double*)exa_malloc(pl->eiLength * sizeof(double));
      
      tr->partitionData[model].substRates        = (double *)exa_malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[model].frequencies       = (double*)exa_malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].empiricalFrequencies       = (double*)exa_malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].tipVector         = (double *)exa_malloc_aligned(pl->tipVectorLength * sizeof(double));
      tr->partitionData[model].symmetryVector    = (int *)exa_malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[model].frequencyGrouping = (int *)exa_malloc(pl->frequencyGroupingLength  * sizeof(int));
      
      tr->partitionData[model].perSiteRates      = (double *)exa_malloc(sizeof(double) * tr->maxCategories);
            
      tr->partitionData[model].nonGTR = FALSE;            

      tr->partitionData[model].gammaRates = (double*)exa_malloc(sizeof(double) * 4);
      tr->partitionData[model].yVector = (unsigned char **)exa_malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

      
      tr->partitionData[model].xVector = (double **)exa_malloc(sizeof(double*) * tr->mxtips);   
      	
      for(j = 0; j < (size_t)tr->mxtips; j++)	        	  	  	  	 
	  tr->partitionData[model].xVector[j]   = (double*)NULL;   

      tr->partitionData[model].xSpaceVector = (size_t *)exa_calloc(tr->mxtips, sizeof(size_t));  

      tr->partitionData[model].sumBuffer = (double *)exa_malloc_aligned(width *
									   (size_t)(tr->partitionData[model].states) *
									   discreteRateCategories(tr->rateHetModel) *
									   sizeof(double));
	    
      tr->partitionData[model].wgt = (int *)exa_malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using exa_calloc() at start up there is only one rate category 0 for all sites */

      tr->partitionData[model].rateCategory = (int *)exa_calloc(width, sizeof(int));

      if(width > 0 && tr->saveMemory)
	{
	  tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	    
	  tr->partitionData[model].gapVector = (unsigned int*)exa_calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));	  	    	  	  
	    
	  tr->partitionData[model].gapColumn = (double *)exa_malloc_aligned(((size_t)tr->mxtips) *								      
									       ((size_t)(tr->partitionData[model].states)) *
									       discreteRateCategories(tr->rateHetModel) * sizeof(double));
	}
      else
	{
	   tr->partitionData[model].gapVectorLength = 0;
	    
	   tr->partitionData[model].gapVector = (unsigned int*)NULL; 	  	    	   
	    
	   tr->partitionData[model].gapColumn = (double*)NULL;	    	    	   
	}              
    }

        
  for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
    myLength += tr->partitionData[model].width;         
   
  /* assign local memory for storing sequence data */

  tr->y_ptr = (unsigned char *)exa_malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);
   
  for(i = 0; i < (size_t)tr->mxtips; i++)
    {
      for( model = 0,  countOffset = 0; model < (size_t)getNumberOfPartitions(tr); model++)
	{
	  tr->partitionData[model].yVector[i+1]   = &tr->y_ptr[i * myLength + countOffset];
	  countOffset +=  tr->partitionData[model].width;
	}
      assert(countOffset == myLength);
    }

  /* figure in data */

  if(tr->manyPartitions)
    {
      for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
	{
	  if(isThisMyPartition(tr, processID, model))
	    {
	      width = tr->partitionData[model].upper - tr->partitionData[model].lower;	     
	      
	      memcpy(&(tr->partitionData[model].wgt[0]), &(tr->aliaswgt[tr->partitionData[model].lower]), sizeof(int) * width);
	    }
	}
    }
  else
    {
      size_t 	   
	globalCounter, 
	r, 
	localCounter;
      
      for(model = 0, globalCounter = 0; model < (size_t)getNumberOfPartitions(tr); model++)
	{
	  for(localCounter = 0, r = (size_t)tr->partitionData[model].lower;  r < (size_t)tr->partitionData[model].upper; r++)
	    {
	      if(r % (size_t)processes == (size_t)processID)
		{
		  tr->partitionData[model].wgt[localCounter] = tr->aliaswgt[globalCounter];	      	     		 		  					     
		  
		  localCounter++;
		}
	      globalCounter++;
	    }
	  assert(localCounter == tr->partitionData[model].width);
	}   
      assert(globalCounter == tr->originalCrunchedLength);
    }
   
  y = (unsigned char *)exa_malloc(sizeof(unsigned char) * tr->originalCrunchedLength);

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength), byteFile);
	
      if(tr->manyPartitions)
	{
	  for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
	    {
	      if(isThisMyPartition(tr, processID, model))	  
		{
		  memcpy(tr->partitionData[model].yVector[i], &(y[tr->partitionData[model].lower]), sizeof(unsigned char) * tr->partitionData[model].width);					    
		  assert(tr->partitionData[model].width == tr->partitionData[model].upper - tr->partitionData[model].lower);
		}
	      else
		assert(tr->partitionData[model].width == 0);
	    }
	}
      else
	{
	  size_t 	  
	    globalCounter, 
	    r, 
	    localCounter;

	  for(model = 0, globalCounter = 0; model < (size_t)getNumberOfPartitions(tr); model++)
	    {
	      for(localCounter = 0, r = (size_t)tr->partitionData[model].lower;  r < (size_t)tr->partitionData[model].upper; r++)
		{
		  if(r % (size_t)processes == (size_t)processID)
		    {		      
		      tr->partitionData[model].yVector[i][localCounter] = y[globalCounter]; 	     
		      
		      localCounter++;
		    }
		  globalCounter++;
		}
	      
	      assert(localCounter == tr->partitionData[model].width);
	    }

	  assert(globalCounter == tr->originalCrunchedLength);
	}
    }

  exa_free(y);
    
  /* initialize gap bit vectors at tips when memory saving option is enabled */
  
  if(tr->saveMemory)
    {
      for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
	{
	  int        
	    undetermined = getUndetermined(tr->partitionData[model].dataType);
	  	 
	  width =  tr->partitionData[model].width;
	    
	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(tr->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(tr->partitionData[model].yVector[j][i] == undetermined)
		    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];	    
	    }     
	}
    }
}



void initializeTree(tree *tr, analdef *adef)
{
  size_t 
    i,
    model;
  
  FILE 
    *byteFile = fopen(byteFileName, "rb");

  double 
    **empiricalFrequencies;	 

  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(size_t), 1, byteFile);
  myBinFread(&(tr->NumberOfModels),         sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);

  empiricalFrequencies = (double **)exa_malloc(sizeof(double *) * getNumberOfPartitions(tr));

  
  if(adef->perGeneBranchLengths)
    tr->numBranches = getNumberOfPartitions(tr);
  else
    tr->numBranches = 1;
  
  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed ExaML-specific */
 
  tr->aliaswgt                   = (int *)exa_malloc(tr->originalCrunchedLength * sizeof(int));
  myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       
  
  tr->rateCategory    = (int *)    exa_calloc(tr->originalCrunchedLength, sizeof(int));	  
  tr->patrat          = (double*)  exa_malloc(tr->originalCrunchedLength * sizeof(double));
  tr->patratStored    = (double*)  exa_malloc(tr->originalCrunchedLength * sizeof(double)); 
  tr->lhs             = (double*)  exa_malloc(tr->originalCrunchedLength * sizeof(double)); 
  
  tr->executeModel   = (boolean *)exa_malloc(sizeof(boolean) * getNumberOfPartitions(tr));
  
  for(i = 0; i < (size_t)getNumberOfPartitions(tr); i++)
    tr->executeModel[i] = TRUE;
   
  setupTree(tr); 


  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int 
	len;
      
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr->nameList[i] = (char*)exa_malloc(sizeof(char) * len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      addword(tr->nameList[i], tr->nameHash, i);        
    }  
 
  for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
    { 
      int 
	len;
      
      pInfo 
	*p = &(tr->partitionData[model]);	   
      
      myBinFread(&(p->states),             sizeof(int), 1, byteFile);
      myBinFread(&(p->maxTipStates),       sizeof(int), 1, byteFile);
      myBinFread(&(p->lower),              sizeof(size_t), 1, byteFile);
      myBinFread(&(p->upper),              sizeof(size_t), 1, byteFile);
      myBinFread(&(p->width),              sizeof(size_t), 1, byteFile);
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
      p->partitionName = (char*)exa_malloc(sizeof(char) * len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);
      
      empiricalFrequencies[model] = (double *)exa_malloc(sizeof(double) * p->states);
      myBinFread(empiricalFrequencies[model], sizeof(double), p->states, byteFile); 
    }     
  
  initializePartitions(tr, byteFile);

  fclose(byteFile);

  initModel(tr, empiricalFrequencies); 
 
  for(model = 0; model < (size_t)getNumberOfPartitions(tr); model++)
    exa_free(empiricalFrequencies[model]);

  exa_free(empiricalFrequencies);
}




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

  if(processID == 0)  
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
	
  if(processID == 0)
    accumulatedTime = 0.0;
	
  /* get the starting tree: here we just parse the tree passed via the command line 
     and do an initial likelihood computation traversal 
     which we maybe should skeip, TODO */
	     
  /* if(tree_file  != NULL) */
  /*   getStartingTree(tr);      */

	   	          
  /* 
     here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
     This should basically be the first call to the library that actually computes something :-)
  */
      
  /* evaluateGeneric(tr, tr->start, TRUE); */
	
  /* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */
      	
  /* treeEvaluate(tr, 1); */

  /* now start the ML search algorithm */
  /* mcmc( tr, adef ); */
  
  exa_main(tr,adef);
  
  
  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
