#include "mem_alloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include "genericParallelization.h"
#include "axml.h"
#include "mem_alloc.h"
/** @file genericParallelization.c
    
    @brief Generic master-worker parallelization with either pthreads or MPI. 
    
    Worker threads/processes mostly work on a local
    tree. Implementationwise, MPI operations are abstracted as good as
    possible via defines (that translate to no-ops or memcpy-calls in
    the pthreads version).

    @todo the code still contains many memory copy operations that
    could be executed more efficiently in-place  
*/



extern unsigned int* mask32; 
extern volatile int jobCycle; 
extern volatile int threadJob; 
extern boolean treeIsInitialized; 

#ifdef MEASURE_TIME_PARALLEL
extern double masterTimePerPhase; 
double timeBuffer[NUM_PAR_JOBS]; 
double timePerRegion[NUM_PAR_JOBS]; 
#endif

extern void initializePartitionData(tree *localTree, partitionList *localPr);
extern void initMemorySavingAndRecom(tree *tr); 
extern char* getJobName(int tmp); 

extern double *globalResult; 
extern volatile char *barrierBuffer;


#ifdef _FINE_GRAIN_MPI
extern MPI_Datatype TRAVERSAL_MPI; 


/** @brief pthreads helper function for adding bytes to buffer.    
 */ 
inline char* addBytes(char *buf, void *toAdd, int numBytes)
{
  memcpy(buf, toAdd, numBytes);  
  return buf + numBytes;  
}

/** @brief pthreads helper function for removing byets from buffer. 
 */ 
inline char* popBytes(char *buf, void *result, int numBytes)
{
  memcpy(result, buf, numBytes); 
  return buf + numBytes;   
}


/** @brief Sets up the MPI environment. 
 *  @param argc   initial argc from main
 *  @param argv   initial argv from main
 */
void initMPI(int argc, char *argv[])
{  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  if(MASTER_P)
    printf("\nThis is RAxML Process Number: %d (MASTER)\n", processID);   
  MPI_Barrier(MPI_COMM_WORLD);
}


/** @brief Traps worker threads.    
    
    @note only relevant for the phtreads version. The master thread
    passes through.

    @param tr the tree 

 */ 
boolean workerTrap(tree *tr, partitionList *pr)
{
  /// @note for the broadcasting, we need to, if the tree structure has already been initialized 
  treeIsInitialized = FALSE; 

  if(NOT MASTER_P) 
    {
      threadData tData; 
      tData.tr = tr; 
      tData.threadNumber = processID;
      tData.pr = pr;
      
      likelihoodThread(&tData);
      return TRUE; 
    }
  return FALSE; 
}


#define ELEMS_IN_TRAV_INFO  9
/** @brief Creates a datastructure for sending the traversal descriptor.
    
    @note This seems to be a very safe method to define your own mpi
   datatypes (often there are problems with padding). But it is not
   entirely for the weak of heart...
 */ 
void defineTraversalInfoMPI()
{
  MPI_Datatype *result  = &TRAVERSAL_MPI; 

  int i ; 
  MPI_Aint base; 
  int blocklen[ELEMS_IN_TRAV_INFO+1] = {1, 1, 1, 1, NUM_BRANCHES, NUM_BRANCHES, 1,1,1,1}; 
  MPI_Aint disp[ELEMS_IN_TRAV_INFO+1];
  MPI_Datatype type[ELEMS_IN_TRAV_INFO+1] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_UB}; 
  traversalInfo desc[2]; 

  MPI_Get_address( desc, disp);
  MPI_Get_address( &(desc[0].pNumber), disp + 1 );
  MPI_Get_address( &(desc[0].qNumber), disp + 2 );  
  MPI_Get_address( &(desc[0].rNumber), disp + 3); 
  MPI_Get_address( desc[0].qz, disp + 4 );
  MPI_Get_address( desc[0].rz, disp + 5 );
  MPI_Get_address( &(desc[0].slot_p), disp + 6);
  MPI_Get_address( &(desc[0].slot_q), disp + 7);
  MPI_Get_address( &(desc[0].slot_r), disp + 8);
  MPI_Get_address( desc + 1, disp + 9);

  base = disp[0]; 
  for(i = 0; i < ELEMS_IN_TRAV_INFO+1; ++i)
    disp[i] -= base;

  MPI_Type_create_struct( ELEMS_IN_TRAV_INFO+1 , blocklen, disp, type, result);
  MPI_Type_commit(result);
}


#endif


/********************/
/* PTHREAD-SPECIFIC */
/********************/
#ifdef _USE_PTHREADS

#ifndef _PORTABLE_PTHREADS
/** @brief Pins a thread to a core (for efficiency). 
    @param tid the thread id
 */ 
void pinToCore(int tid)
{
  cpu_set_t cpuset;

  CPU_ZERO(&cpuset);    
  CPU_SET(tid, &cpuset);

  if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
    {
      printBothOpen("\n\nThere was a problem finding a physical core for thread number %d to run on.\n", tid);
      printBothOpen("Probably this happend because you are trying to run more threads than you have cores available,\n");
      printBothOpen("which is a thing you should never ever do again, good bye .... \n\n");
      assert(0);
    }
}
#endif

/**  @brief Starts the worker threads. 
     @param tr
 */ 
void startPthreads(tree *tr, partitionList *pr)
{
  pthread_t *threads;
  pthread_attr_t attr;
  int rc, t;
  threadData *tData;
  treeIsInitialized = FALSE; 

  jobCycle        = 0;
  threadJob       = 0;

  printf("\nThis is the RAxML Master Pthread\n");  

#if (NOT defined(_USE_PTHREADS) && defined( MEASURE_TIME_PARALLEL))
  timeBuffer = rax_calloc(NUM_PAR_JOBS * tr->numberOfThreads, sizeof(double)); 
#endif

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  threads    = (pthread_t *)rax_malloc((size_t)tr->numberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)rax_malloc((size_t)tr->numberOfThreads * sizeof(threadData));

  barrierBuffer            = (volatile char *)  rax_malloc(sizeof(volatile char)   *  (size_t)tr->numberOfThreads);

  for(t = 0; t < tr->numberOfThreads; t++)
    barrierBuffer[t] = 0;

  for(t = 1; t < tr->numberOfThreads; t++)
    {
      tData[t].tr  = tr;
      tData[t].pr  = pr;
      tData[t].threadNumber = t;
      rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t]));
      if(rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
    }
}
#endif

#ifdef MEASURE_TIME_PARALLEL
static void reduceTimesWorkerRegions(tree *tr, double *mins, double *maxs)
{
  int tid = tr->threadID; 
  int i,j ; 
  double reduction[NUM_PAR_JOBS * tr->numberOfThreads]; 

  ASSIGN_GATHER(reduction, timeBuffer, NUM_PAR_JOBS, DOUBLE, tr->threadID); 

#ifdef _USE_PTHREADS
  /* we'd need a proper barrier here... this evaluation is mostly interesting for MPI  */
  printf("\n\ncomment out MEASURE_TIME_PARALLEL\n\n");   
  assert(0); 
#else 
   MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* find min and max time */
  if(MASTER_P)
    {      
      for(j = 0; j < NUM_PAR_JOBS; ++j)
	{
	  boolean isFirst = TRUE; 
	  for(i = 0; i < tr->numberOfThreads; ++i)
	    {	      
	      double num = timeBuffer[i * NUM_PAR_JOBS + j]; 
	      if(isFirst || num < mins[j])
		mins[j] = num; 
	      if(isFirst || num > maxs[j])
		maxs[j] = num; 
	      isFirst = FALSE; 
	    }
	}	
    }  
}


static void printParallelTimePerRegion(double *mins, double *maxs)
{
  int i; 
  double allTime = 0; 
  double relTime[NUM_PAR_JOBS+1]; 
  for(i = 0; i < NUM_PAR_JOBS+1; ++i)
    allTime += timePerRegion[i]; 
  for(i = 0; i < NUM_PAR_JOBS+1; ++i)
    relTime[i] = (timePerRegion[i] / allTime) * 100 ; 

  printf("\n\nTime spent per region \nmasterTimeAbs\tmasterTimeRel\tloadBalance\tcommOverhead\n"); 
  for(i = 0; i < NUM_PAR_JOBS; ++i )
    if(timePerRegion[i] != 0)
      printf("%f\t%.2f\%\t%.2f\%\t%.2f\%\t%s\n", timePerRegion[i], relTime[i], (maxs[i] - mins[i]) * 100  / maxs[i] , (maxs[i] - timePerRegion[i]) * 100  / timePerRegion[i], getJobName(i) ); 
  printf("================\n%f\t%.2f\%\tSEQUENTIAL\n", timePerRegion[NUM_PAR_JOBS], relTime[i]); 
  printf("loadbalance: (minT - maxT) / maxT, \twhere minT is the time the fastest worker took for this region (maxT analogous) \n"); 
  printf("commOverhead: (maxWorker - masterTime) / masterTime, \t where maxWorker is the time the slowest worker spent in this region and masterTime is the time the master spent in this region (e.g., doing additional reduction stuff)\n"); 
}
#endif




/* function that computes per-site log likelihoods in pthreads */

/** @brief worker threads evaluate the likelihood on their sites. 
    @param tr the tree
    @param likelihood array (?)  
    @param n number of threads
    @param tid thread id 
 */ 
void perSiteLogLikelihoodsPthreads(tree *tr, partitionList *pr, double *lhs, int n, int tid)
{
  size_t 
    model, 
    i;

  for(model = 0; model < pr->numberOfPartitions; model++)
    {      
      size_t 
	localIndex = 0;

      /* decide if this partition is handled by the thread when -Q is ativated 
	 or when -Q is not activated figure out which sites have been assigned to the 
	 current thread */

      boolean 
	execute = ((tr->manyPartitions && isThisMyPartition(pr, tid, model)) || (!tr->manyPartitions));

      /* if the entire partition has been assigned to this thread (-Q) or if -Q is not activated 
	 we need to compute some per-site log likelihoods with thread tid for this partition */

      if(execute)
	for(i = pr->partitionData[model]->lower;  i < pr->partitionData[model]->upper; i++)
	  {
	    /* if -Q is active we compute all per-site log likelihoods for the partition,
	       othwerise we only compute those that have been assigned to thread tid 
	       using the cyclic distribution scheme */

	    if(tr->manyPartitions || (i % n == tid))
	      {
		double 
		  l;

		/* now compute the per-site log likelihood at the current site */

		switch(tr->rateHetModel)
		  {
		  case CAT:
		    l = evaluatePartialGeneric (tr, pr, localIndex, pr->partitionData[model]->perSiteRates[pr->partitionData[model]->rateCategory[localIndex]], model);
		    break;
		  case GAMMA:
		    l = evaluatePartialGeneric (tr, pr, localIndex, 1.0, model);
		    break;
		  default:
		    assert(0);
		  }

		/* store it in an array that is local in memory to the current thread,
		   see function collectDouble() in axml.c for understanding how we then collect these 
		   values stored in local arrays from the threads */

		lhs[i] = l;

		localIndex++;
	      }
	  }
    }
}

/** @brief Check, if partition is assign to this worker.
    @param localTree the local tree 
    @param tid the thread id
    @param model the partition id
 */ 
boolean isThisMyPartition(partitionList *localPr, int tid, int model)
{ 
  if(localPr->partitionData[model]->partitionAssignment == tid)
    return TRUE;
  else
    return FALSE;
}

/** @brief Computes partition size for all partitions (in case full partitions are assigns to workers). 

    @param localTree the local tree 
    
    @param tid thread id    
 */ 
void computeFractionMany(partitionList *localPr, int tid)
{
  int
    sites = 0;

  int   
    model;

  for(model = 0; model < localPr->numberOfPartitions; model++)
    {
      if(isThisMyPartition(localPr, tid, model))
	{	 
    	  localPr->partitionData[model]->width = localPr->partitionData[model]->upper - localPr->partitionData[model]->lower;
	  sites += localPr->partitionData[model]->width;
	}
      else       	  
    	  localPr->partitionData[model]->width = 0;
    }


}


/** @brief Computes partition size for all partitions (for cyclic distribution of sites)
    
    @param localTree the local tree 
    @param tid thread id
    @param n number of workers
 */ 
void computeFraction(partitionList *localPr, int tid, int n)
{
  int
    i,
    model;

  for(model = 0; model < localPr->numberOfPartitions; model++)
    {
      int width = 0;

      for(i = localPr->partitionData[model]->lower; i < localPr->partitionData[model]->upper; i++)
	if(i % n == tid)
	  width++;
      localPr->partitionData[model]->width = width;
    }
}



/** @brief Compare partition sizes. 
    @param p1 pointer to a partition
    @param p2 pointer to another partition
 */ 
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


/** @brief Top-level function for the multi processor scheduling scheme (assigns full partitions to workers). 
    
   tr->manyPartitions is set to TRUE if the user has indicated via -Q
   that there are substantially more partitions than threads/cores
   available. In that case we do not distribute sites from each
   partition in a cyclic fashion to the cores , but distribute entire
   partitions to cores.  Achieving a good balance of alignment sites
   to cores boils down to the multi-processor scheduling problem known
   from theoretical comp. sci.  which is NP-complete.  We have
   implemented very simple "standard" heuristics for solving the
   multiprocessor scheduling problem that turn out to work very well
   and are cheap to compute.
   
   @param tr tree 
   @param worker id 
*/

void multiprocessorScheduling(tree *tr, partitionList *pr, int tid)
{
  int 
    s,
    model,
    modelStates[2] = {4, 20},
    numberOfPartitions[2] = {0 , 0},
      arrayLength = sizeof(modelStates) / sizeof(int);

      /* check that we have not addedd any new models for data types with a different number of states
	 and forgot to update modelStates */

      // TODO: Check whether we need this line!
      // tr->partitionAssignment = (int *)rax_malloc((size_t)tr->NumberOfModels * sizeof(int));

      for(model = 0; model < pr->numberOfPartitions; model++)
	{        
	  boolean 
	    exists = FALSE;

	  for(s = 0; s < arrayLength; s++)
	    {
	      exists = exists || (pr->partitionData[model]->states == modelStates[s]);
	      if(pr->partitionData[model]->states == modelStates[s])
		numberOfPartitions[s] += 1;
	    }

	  assert(exists);
	}

      if(tid == 0)
	printBothOpen("\nMulti-processor partition data distribution enabled (\"-Q\" option)\n");

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
#ifndef _FINE_GRAIN_MPI
		n = tr->numberOfThreads,
#else
		n = processes,
#endif
		p = numberOfPartitions[s],    
		*assignments = (int *)rax_calloc((size_t)n, sizeof(int));  

	      partitionType 
		*pt = (partitionType *)rax_malloc(sizeof(partitionType) * (size_t)p);



	      for(i = 0, k = 0; i < pr->numberOfPartitions; i++)
		{
		  if(pr->partitionData[i]->states == modelStates[s])
		    {
		      pt[k].partitionNumber = i;
		      pt[k].partitionLength = pr->partitionData[i]->upper - pr->partitionData[i]->lower;
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
		  assert(pt[i].partitionNumber >= 0 && pt[i].partitionNumber < pr->numberOfPartitions);
		  pr->partitionData[pt[i].partitionNumber]->partitionAssignment = minIndex;
		}

	      if(tid == 0)
		{
		  for(i = 0; i < n; i++)	       
		    printBothOpen("Process %d has %d sites for %d state model \n", i, assignments[i], modelStates[s]);		  		

		  printBothOpen("\n");
		}

	      for(i = 0; i < n; i++)
		checkSum += (size_t)assignments[i];

	      assert(sum == checkSum);

	      rax_free(assignments);
	      rax_free(pt);
	    }
	}
}



/** @brief Reduce the first and second derivative of the likelihood
    function.
    
    We collect the first and second derivatives from the various
    threads and sum them up. It's similar to what we do in
    evaluateGeneric() with the only difference that we have to collect
    two values (firsrt and second derivative) instead of onyly one (the
    log likelihood

   @warning operates on global reduction buffers \a globalResult
   
   @param tr tree 
   @param dlnLdlz first derivative
   @param d2lnLdlz2 second derivative
*/
void branchLength_parallelReduce(tree *tr, double *dlnLdlz,  double *d2lnLdlz2, int numBranches )
{
#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS

  /* only the master executes this  */
  assert(tr->threadID == 0); 
  
  int b; 
  int t; 
  for(b = 0; b < numBranches; ++b)
    {
      dlnLdlz[b] = 0; 
      d2lnLdlz2[b] = 0; 

      for(t = 0; t < tr->numberOfThreads; ++t)
	{
	  dlnLdlz[b] += globalResult[t * numBranches * 2 + b ];
	  d2lnLdlz2[b] += globalResult[t * numBranches * 2 + numBranches + b];
	}
    }
#else 
  memcpy(dlnLdlz, globalResult, sizeof(double) * numBranches);
  memcpy(d2lnLdlz2, globalResult + numBranches, sizeof(double) * numBranches);
#endif
}



/** @brief Read from buffer or writes rates into buffer.  Return
    number of elems written.
    
   @param buf the buffer
   @param srcTar pointer to either source or destination array  
   @param tr tree 
   @param n number of workers 
   @param tid process id 
   @param read TRUE, if read-mode  
   @param countOnly  if TRUE, simply return the number of elements 

*/
static int doublesToBuffer(double *buf, double *srcTar, tree *tr, partitionList *pr, int n, int tid, boolean read, boolean countOnly)
{
  int 
    model,
    i;
  double 
    *initPtr = buf; 

  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      if(tr->manyPartitions)
	{
	  if(isThisMyPartition(pr, tid, model))
	    for(i = pr->partitionData[model]->lower; i < pr->partitionData[model]->upper; i++)
	      {
		if(NOT countOnly)
		  {
		    if(read)
		      *buf = srcTar[i]; 
		    else 
		      srcTar[i] = *buf; 
		  }
		buf++;
	      }	  
	}
      
      else
	{
	  for(i = pr->partitionData[model]->lower; i < pr->partitionData[model]->upper; i++)
	    if(i % n == tid)
	      {
		if(NOT countOnly)
		  {
		    if(read)
		      *buf = srcTar[i];
		    else 
		      srcTar[i] = *buf; 
		  }
		buf++; 
	      }
	}
    }
  
  return buf - initPtr; 
}


/** @brief broadcast rates after rate optimization. 
    
    @param tr tree 
    @param localTree local tree 
    @param n number of workers 
    @param tid worker id 
    
    @todo mpi_alltoallv/w may be more efficient, but it is a hell to set up
 */ 
void broadcastAfterRateOpt(tree *tr, tree *localTree, partitionList *pr, int n, int tid)
{				  
  int
    num1 = 0,
    num2 = 0,
    num3 = 0, 
    i ; 
    
  for(i = 0; i < n; ++i)
    {
      double
	allBuf[tr->originalCrunchedLength * 3],
	buf1[tr->originalCrunchedLength],
	buf2[tr->originalCrunchedLength], 
	buf3[tr->originalCrunchedLength]; 

#ifdef _USE_PTHREADS
      if(i != tid)
	continue; 
#endif
      int numDouble = 0; 
      
      /* extract doubles  */

      num1 = doublesToBuffer(buf1, localTree->patrat, tr, pr, n,i, TRUE, i!= tid);
      num2 = doublesToBuffer(buf2, localTree->patratStored, tr, pr, n,i, TRUE, i!= tid);
      num3 = doublesToBuffer(buf3, localTree->lhs, tr, pr, n,i, TRUE, i!= tid);

      numDouble += num1 + num2 + num3; 

      /* copy doubles  */
      
      memcpy(allBuf, buf1, num1 * sizeof(double)); 
      memcpy(allBuf + num1, buf2, num2 * sizeof(double)); 
      memcpy(allBuf + (num1 + num2) , buf3, num3 * sizeof(double)); 

      BCAST_BUF(allBuf, numDouble, MPI_DOUBLE, i); 

      memcpy(buf1, allBuf, num1 * sizeof(double)); 
      memcpy(buf2, allBuf + num1, num2 * sizeof(double)); 
      memcpy(buf3, allBuf + (num1 + num2), num3 * sizeof(double)); 
      
      /* re-insert doubles  */
      int assertCtr = 0; 
      assertCtr += doublesToBuffer(buf1, tr->patrat, tr, pr, n,i,FALSE, FALSE);
      assertCtr += doublesToBuffer(buf2, tr->patratStored, tr, pr, n,i,FALSE, FALSE);
      assertCtr += doublesToBuffer(buf3, tr->lhs, tr, pr, n,i,FALSE, FALSE);

      assert(assertCtr == numDouble); 
    }
}


/** @brief Collect doubles from workers to master.
    
    @param dst destination array
    @param src source array
    @param tr tree 
    @param n number of workers 
    @param tid worker id 
 */
static void collectDouble(double *dst, double *src, tree *tr, partitionList *pr, int n, int tid)
{
  int i; 
  double 
    resultBuf[tr->originalCrunchedLength],
    buf[tr->originalCrunchedLength]; 
  int
    assertNum, 
    displacements[tr->numberOfThreads]; 


  /* gather own values into buffer  */
  int numberCollected = doublesToBuffer(buf, src, tr, pr,n,tid,TRUE, FALSE);

#ifdef _FINE_GRAIN_MPI 
  /* this communicates all the values to the master */
  
  int numberPerWorker[tr->numberOfThreads];     
  if(MASTER_P)			/* master counts number to receive, receives and writes back */
    {
      for(i = 0; i < n; ++i)
	{
	  numberPerWorker[i] = doublesToBuffer(buf,src,tr,pr,n,i,FALSE, TRUE);
	  displacements[i] = i == 0 ? 0 : displacements[i-1] + numberPerWorker[i-1]; 
	}
      
      MPI_Gatherv(buf, numberCollected, MPI_DOUBLE,
		  resultBuf, numberPerWorker, displacements,  MPI_DOUBLE,
		  0, MPI_COMM_WORLD); 

      double *bufPtr = resultBuf; 
      for(i = 0 ; i < n; ++i)
	{
	  int numberWritten = doublesToBuffer(bufPtr, dst,tr,pr,n,i, FALSE, FALSE);
	  bufPtr += numberWritten; 
	  assertNum += numberWritten; 
	}    
      
      assert(assertNum == tr->originalCrunchedLength);
    }
  else 				/* workers only send their buffer   */
    MPI_Gatherv(buf, numberCollected, MPI_DOUBLE, resultBuf, numberPerWorker, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
#else 
  /* pthread only writes to global space  */  
  assertNum = doublesToBuffer(buf, dst,tr,pr,n,tid, FALSE, FALSE);
  assert(assertNum == numberCollected); 
#endif
}

/** @brief broadcast a new alpha (for the GAMMA model)
    @param localTree local tree 
    @param tr tree 
    @param tid worker id 
 */
static void broadCastAlpha(tree *localTree, tree *tr, partitionList *localPr, partitionList *pr, int tid)
{
  int 
    model; 

  int
    i,
    bufSize = localPr->numberOfPartitions * 4 * sizeof(double);

  char bufDbl[bufSize], 
    *bufPtrDbl = bufDbl;   

  RECV_BUF(bufDbl, bufSize, MPI_BYTE); 

  for(model = 0; model < localPr->numberOfPartitions; model++)
    for(i = 0; i < 4; ++i)
      ASSIGN_BUF_DBL(localPr->partitionData[model]->gammaRates[i], pr->partitionData[model]->gammaRates[i]);
  
  SEND_BUF(bufDbl, bufSize, MPI_BYTE);  
}


/** @brief Master broadcasts rates.
    
    @param localTree local tree 
    @param tr tree 
    @param tid worker id     
 */ 
static void broadCastRates(tree *localTree, tree *tr, partitionList *localPr, partitionList *pr, int tid)
{
  int 
    model;

  /* determine size of buffer needed first */
  int bufSize = 0; 
#ifdef _FINE_GRAIN_MPI
  for(model = 0; model < localPr->numberOfPartitions; ++model )
    {	  
      const partitionLengths *pl = getPartitionLengths(pr->partitionData[model]); /* this is constant, isnt it?  */
      bufSize += (pl->eignLength + pl->evLength + pl->eiLength + pl->tipVectorLength) * sizeof(double) ; 
    }
#endif      
 
  char bufDbl[bufSize], 
    *bufPtrDbl = bufDbl; 
  RECV_BUF(bufDbl, bufSize, MPI_BYTE);
  int i ; 

  for(model = 0; model < localPr->numberOfPartitions; model++)
    {
      const partitionLengths *pl = getPartitionLengths(pr->partitionData[model]); /* this is constant, isnt it?  */

      for(i = 0; i < pl->eignLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->EIGN[i], pr->partitionData[model]->EIGN[i]);
      for(i = 0; i < pl->evLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->EV[i],pr->partitionData[model]->EV[i]);
      for(i = 0; i  < pl->eiLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->EI[i], pr->partitionData[model]->EI[i]);
      for(i = 0; i < pl->tipVectorLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->tipVector[i],   pr->partitionData[model]->tipVector[i]);
    }
  SEND_BUF(bufDbl, bufSize, MPI_BYTE); /*  */
}


/** @brief likelihood evaluation call with subsequent reduce operation. 

    @param localTree local tree 
    @param tid worker id 
 */ 
static void reduceEvaluateIterative(tree *localTree, partitionList *localPr, int tid)
{
  int model;

  evaluateIterative(localTree, localPr);

  /* when this is done we need to write the per-thread log likelihood to the 
     global reduction buffer. Tid is the thread ID, hence thread 0 will write its 
     results to reductionBuffer[0] thread 1 to reductionBuffer[1] etc.

     the actual sum over the entries in the reduction buffer will then be computed 
     by the master thread which ensures that the sum is determinsitic */

#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS
  /* 
     aberer: I implemented this as a mpi_gather operation into this buffer, 
     pthreads version emulates this gather; 
     master takes care of the reduction; 
  */

  double buf[localPr->numberOfPartitions];
  for(model = 0; model < localPr->numberOfPartitions; ++model)
    buf[model] = localPr->partitionData[model]->partitionLH;

  /* either make reproducible or efficient */
  ASSIGN_GATHER(globalResult, buf, localPr->numberOfPartitions, DOUBLE, tid);
#else 
  /* the efficient mpi version: a proper reduce  */
  double buf[localPr->numberOfPartitions];
  for(model = 0; model < localPr->numberOfPartitions; ++model)
      buf[model] = localPr->partitionData[model]->partitionLH;
  double targetBuf[localPr->numberOfPartitions];
  memset(targetBuf, 0, sizeof(double) * localPr->numberOfPartitions);
  MPI_Reduce(buf, targetBuf, localPr->numberOfPartitions, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(MASTER_P) {
	  for(model = 0; model < localPr->numberOfPartitions; ++model) {
		  localPr->partitionData[model]->partitionLH = targetBuf[model];
	  }
  }
#endif
}



/*@ @brief Broadcast the traversal descriptor to worker threads. 

  The one below is a hack we are re-assigning the local pointer to
  the global one the memcpy version below is just for testing and
  preparing the fine-grained MPI BlueGene version

  @param localTree local tree 
  @param tr tree 
*/
/* TODO: we should reset this at some point, the excplicit copy is just done for testing */
inline static void broadcastTraversalInfo(tree *localTree, tree *tr, partitionList *localPr, partitionList *pr)
{
  /* @todo these two regions could be joined */
#ifdef _USE_PTHREADS
  /* memcpy -> memmove (see ticket #43). This function is sometimes called with localTree == tr,
   * in which case some memcpy implementations can corrupt the buffers.
   */
  
  localTree->td[0].functionType =            tr->td[0].functionType;
  localTree->td[0].count =                   tr->td[0].count ;
  localTree->td[0].traversalHasChanged =     tr->td[0].traversalHasChanged;

  memmove(localTree->td[0].executeModel,    tr->td[0].executeModel,    sizeof(boolean) * localPr->numberOfPartitions);
  memmove(localTree->td[0].parameterValues, tr->td[0].parameterValues, sizeof(double) * localPr->numberOfPartitions);
  
  if(localTree->td[0].traversalHasChanged)
    memmove(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));

#else
  /* MPI */
  /* like in raxml-light: first we send a small message, if the
     travesalDescriptor is longer, then resend */
  
  int length = treeIsInitialized ? localPr->numberOfPartitions : 0;
  char broadCastBuffer[messageSize(length)]; 
  char *bufPtr = broadCastBuffer; 
  int i; 

  RECV_BUF(broadCastBuffer, messageSize(length), MPI_BYTE); 

  ASSIGN_BUF(localTree->td[0].functionType, tr->td[0].functionType , int);   
  ASSIGN_BUF(localTree->td[0].count,  tr->td[0].count , int); 
  ASSIGN_BUF(localTree->td[0].traversalHasChanged, tr->td[0].traversalHasChanged , int); 

  if(treeIsInitialized)  
    { 
      for(i = 0; i < localPr->numberOfPartitions; ++i)
	{
	  ASSIGN_BUF(localTree->td[0].executeModel[i],      tr->td[0].executeModel[i], int); 
	  ASSIGN_BUF(localTree->td[0].parameterValues[i],	 tr->td[0].parameterValues[i], double); 
	}      

      for(i = 0; i < TRAVERSAL_LENGTH; ++i )
	ASSIGN_BUF(localTree->td[0].ti[i], tr->td[0].ti[i], traversalInfo); 
    }
    
  SEND_BUF(broadCastBuffer, messageSize(length), MPI_BYTE); 

  /* now we send the second part of the traversal descriptor, if we
     exceed the pre-set number of elements */
  if(treeIsInitialized && localTree->td[0].count > TRAVERSAL_LENGTH) 
    {
      /* lets use the MPI_Datatype for this thing, what I've read it's
	 supposed to be more secure and efficient */
      MPI_Bcast(localTree->td[0].ti + TRAVERSAL_LENGTH, localTree->td[0].count - TRAVERSAL_LENGTH, TRAVERSAL_MPI, 0, MPI_COMM_WORLD );
    }
#endif
}


/** @brief helper that yields a string representation of a parallel region. 
    
    @param type type of parallel region
 */ 
char* getJobName(int type)
{
  switch(type)  
    {
    case  THREAD_NEWVIEW:       
      return "THREAD_NEWVIEW";
    case THREAD_EVALUATE: 
      return "THREAD_EVALUATE";
    case THREAD_MAKENEWZ: 
      return "THREAD_MAKENEWZ";
    case THREAD_MAKENEWZ_FIRST: 
      return "THREAD_MAKENEWZ_FIRST";
    case THREAD_RATE_CATS: 
      return "THREAD_RATE_CATS";
    case THREAD_COPY_RATE_CATS: 
      return "THREAD_COPY_RATE_CATS";
    case THREAD_COPY_INIT_MODEL: 
      return "THREAD_COPY_INIT_MODEL";
    case THREAD_INIT_PARTITION: 
      return "THREAD_INIT_PARTITION";
    case THREAD_OPT_ALPHA: 
      return "THREAD_OPT_ALPHA";
    case THREAD_OPT_RATE: 
      return "THREAD_OPT_RATE";
    case THREAD_COPY_ALPHA: 
      return "THREAD_COPY_ALPHA";
    case THREAD_COPY_RATES: 
      return "THREAD_COPY_RATES";
    case THREAD_PER_SITE_LIKELIHOODS: 
      return "THREAD_PER_SITE_LIKELIHOODS";
    case THREAD_NEWVIEW_ANCESTRAL: 
      return "THREAD_NEWVIEW_ANCESTRAL";
    case THREAD_GATHER_ANCESTRAL: 
      return "THREAD_GATHER_ANCESTRAL";
    case THREAD_EXIT_GRACEFULLY: 
      return "THREAD_EXIT_GRACEFULLY";
    default: assert(0); 
    }
}


/**
   @brief Generic entry point for parallel regions (mostly broadcasts
   traversal descriptor first).

   This function here handles all parallel regions in the Pthreads
   version, when we enter this function masterBarrier() has ben called
   by the master thread from within the sequential part of the
   program, tr is the tree at the master thread, localTree the tree at
   the worker threads

   While this is not necessary, adress spaces of threads are indeed
   separated for easier transition to a distributed memory paradigm
   
   @param tr tree 
   @param localTree local tree 
   @param tid worker id 
   @param n number of workers 
*/
boolean execFunction(tree *tr, tree *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
{
  int
    i,
    model,
    localCounter;

#ifdef MEASURE_TIME_PARALLEL
  double timeForParallelRegion = gettime();
#endif


#ifdef _USE_PTHREADS
  /* some stuff associated with the barrier implementation using Pthreads and busy wait */
  int currentJob = threadJob >> 16;
#endif

  /* here the master sends and all threads/processes receive the traversal descriptor */
  broadcastTraversalInfo(localTree, tr, localPr, pr);

#ifdef _USE_PTHREADS
  /* make sure that nothing is going wrong */
  assert(currentJob == localTree->td[0].functionType);
#else   
  localTree = tr; 
  int currentJob = localTree->td[0].functionType; 
#endif
#ifdef DEBUG_PARALLEL
  printf("[%d] working on %s\n", tid, getJobName(currentJob)); 
#endif  

  switch(currentJob)
    { 
    case THREAD_NEWVIEW: 
      /* just a newview on the fraction of sites that have been assigned to this thread */

      newviewIterative(localTree, localPr, 0);
      break;     
    case THREAD_EVALUATE: 
      reduceEvaluateIterative(localTree, localPr, tid);
      break;	
    case THREAD_MAKENEWZ_FIRST:

      /* this is the first call from within makenewz that requires getting the likelihood vectors to the left and 
         right of the branch via newview and doing some precomputations.
	 
         For details see comments in makenewzGenericSpecial.c 
      */
    case  THREAD_MAKENEWZ:
      {	
	double
	  dlnLdlz[NUM_BRANCHES],
	  d2lnLdlz2[NUM_BRANCHES]; 

	if(localTree->td[0].functionType == THREAD_MAKENEWZ_FIRST)
	  makenewzIterative(localTree, localPr);
	execCore(localTree, localPr, dlnLdlz, d2lnLdlz2);

	/* gather the first and second derivatives that have been written by each thread */
	/* as for evaluate above, the final sum over the derivatives will be computed by the 
	   master thread in its sequential part of the code */

	int numBranches = localPr->perGeneBranchLengths?localPr->numberOfPartitions:1;

#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS
	/* MPI: implemented as a gather again, pthreads: just buffer copying */	
	double buf[ 2 * numBranches];
	memcpy( buf, dlnLdlz, numBranches * sizeof(double) );
	memcpy(buf + numBranches, d2lnLdlz2, numBranches * sizeof(double));

	ASSIGN_GATHER(globalResult, buf,  2 * numBranches, DOUBLE, tid);
#else 	
	double result[numBranches];
	memset(result,0, numBranches * sizeof(double));
	MPI_Reduce( dlnLdlz , result , numBranches, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MASTER_P)
	  memcpy(globalResult, result, sizeof(double) * numBranches);
	
	memset(result,0,numBranches * sizeof(double));
	MPI_Reduce( d2lnLdlz2 , result , numBranches, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MASTER_P)
	  memcpy(globalResult + numBranches, result, sizeof(double) * numBranches);
#endif
      }

      break;

    case THREAD_INIT_PARTITION:       

      /* broadcast data and initialize and allocate arrays in partitions */
      
      initializePartitions(tr, localTree, pr, localPr, tid, n);

      break;          
    case THREAD_COPY_ALPHA: 
    case THREAD_OPT_ALPHA:
      /* this is when we have changed the alpha parameter, inducing a change in the discrete gamma rate categories.
	 this is called when we are optimizing or sampling (in the Bayesioan case) alpha parameter values */
      
      /* distribute the new discrete gamma rates to the threads */
      broadCastAlpha(localTree,tr, localPr,pr, tid);

      /* compute the likelihood, note that this is always a full tree traversal ! */
      if(localTree->td[0].functionType == THREAD_OPT_ALPHA)
	reduceEvaluateIterative(localTree, localPr, tid);

      break;           
    case THREAD_OPT_RATE:
    case THREAD_COPY_RATES:

      /* if we are optimizing the rates in the transition matrix Q this induces recomputing the eigenvector eigenvalue 
	 decomposition and the tipVector as well because of the special numerics in RAxML, the matrix of eigenvectors 
	 is "rotated" into the tip lookup table.

	 Hence if the sequantial part of the program that steers the Q matrix rate optimization has changed a rate we 
	 need to broadcast all eigenvectors, eigenvalues etc to each thread 
      */

      broadCastRates(localTree, tr, localPr, pr, tid);

      /* now evaluate the likelihood of the new Q matrix, this always requires a full tree traversal because the changes need
	 to be propagated throughout the entire tree */

      if(localTree->td[0].functionType == THREAD_OPT_RATE)
	reduceEvaluateIterative(localTree, localPr, tid);

      break;                       
    case THREAD_COPY_INIT_MODEL:
      {

	/* need to be very careful here ! THREAD_COPY_INIT_MODEL is also used when the program is restarted 
	   it is hence not sufficient to just initialize everything by the default values ! */

	broadCastRates(localTree, tr, localPr, pr, tid);
	broadCastAlpha(localTree, tr, localPr, pr, tid); /* isnt that only executed when we are on gamma?  */

	/*
	  copy initial model parameters, the Q matrix and alpha are initially, when we start our likelihood search 
	  set to default values. 
	  Hence we need to copy all those values that are required for computing the likelihood 
	  with newview(), evaluate() and makenez() to the private memory of the threads 
	*/


	if( localTree->rateHetModel == CAT) /* TRICKY originally this should only be executed by workers  */
	  { 	    
	    int bufSize = 2 * localTree->originalCrunchedLength * sizeof(double); 
	    char bufDbl[bufSize], 
	      *bufPtrDbl = bufDbl; 

	    RECV_BUF(bufDbl, bufSize,MPI_BYTE); 

	    /* this should be local  */
	    for(model = 0; model < localPr->numberOfPartitions; model++)
	      localPr->partitionData[model]->numberOfCategories      = pr->partitionData[model]->numberOfCategories;


	    /* this is only relevant for the PSR model, we can worry about this later */
	    for(i = 0; i < localTree->originalCrunchedLength; ++i)
	      {
		ASSIGN_BUF_DBL(localTree->patrat[i], tr->patrat[i]);
		ASSIGN_BUF_DBL(localTree->patratStored[i], tr->patratStored[i]); 
	      }

	    SEND_BUF(bufDbl, bufSize, MPI_BYTE); 
	  }
      } 
      break;    
    case THREAD_RATE_CATS: 
      {
	/* this is for optimizing per-site rate categories under PSR, let's worry about this later */

	ASSIGN_DBL( localTree->lower_spacing,  tr->lower_spacing);
	ASSIGN_DBL( localTree->upper_spacing,  tr->upper_spacing);

	optRateCatPthreads(localTree, localPr, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);

	broadcastAfterRateOpt(tr, localTree, pr, n,  tid);
      }
      break;
    case THREAD_COPY_RATE_CATS:
      {
	/* 
	   this is invoked when we have changed the per-site rate category assignment
	   In essence it distributes the new per site rates to all threads 

	   The pthread-version here simply assigns everything as ought to
	   be. The MPI-version is configured to write to a buffer instead
	   and SEND (master) or RECV (workers) it.

	*/

	/* 
	   start of communication part 
	*/

	int i, 
	  /* buf[localPr->numberOfPartitions], */
	  /* assertCtr = 0,  */
	  dblBufSize = 0; 
	int bufSize = localPr->numberOfPartitions * sizeof(int);
	char buf[bufSize], 
	  *bufPtr = buf; 
     
	RECV_BUF(buf, bufSize, MPI_BYTE);

	for( model = 0; model < localPr->numberOfPartitions; ++model)
	  {
	    ASSIGN_BUF(localPr->partitionData[model]->numberOfCategories, pr->partitionData[model]->numberOfCategories, int);
	    dblBufSize += localPr->partitionData[model]->numberOfCategories * sizeof(double);
	  }

	SEND_BUF(buf, bufSize, MPI_BYTE); 


	dblBufSize += 2 * localTree->originalCrunchedLength * sizeof(double); 

	char bufDbl[dblBufSize],
	  *bufPtrDbl = bufDbl;

	RECV_BUF(bufDbl, dblBufSize, MPI_BYTE); 

	for(i = 0; i < localTree->originalCrunchedLength; ++i)
	  {	 
	    ASSIGN_BUF_DBL(localTree->patrat[i], tr->patrat[i]); 
	    ASSIGN_BUF_DBL(localTree->patratStored[i], tr->patratStored[i]); 
	  }

	for( model = 0; model < localPr->numberOfPartitions; ++model)
	  for(i = 0; i < localPr->partitionData[model]->numberOfCategories; i++)
	    ASSIGN_BUF_DBL(localPr->partitionData[model]->perSiteRates[i], pr->partitionData[model]->perSiteRates[i]);

	SEND_BUF(bufDbl, dblBufSize, MPI_BYTE); 


	/* lets test, if it is a good idea to send around the basic categories  */
#ifdef _FINE_GRAIN_MPI
	/* TODO this is inefficient, but is seems to have a small impact on performance */
	MPI_Bcast(tr->rateCategory, tr->originalCrunchedLength, MPI_INT, 0, MPI_COMM_WORLD); 
#endif


	/* 
	   now re-assign values 
	*/
	for(model = 0; model < localPr->numberOfPartitions; model++)
	  {
	    if(localTree->manyPartitions)
	      {
		if(isThisMyPartition(localPr, tid, model))
		  for(localCounter = 0, i = localPr->partitionData[model]->lower;  i < localPr->partitionData[model]->upper; i++, localCounter++)
		    {	     
		      localPr->partitionData[model]->rateCategory[localCounter] = tr->rateCategory[i];
		    } 
	      }
	    else	  
	      {
		for(localCounter = 0, i = localPr->partitionData[model]->lower;  i < localPr->partitionData[model]->upper; i++)
		  {
		    if(i % n == tid)
		      {		 
			localPr->partitionData[model]->rateCategory[localCounter] = tr->rateCategory[i];

			localCounter++;
		      }
		  }
	      }
	  }
      }
      break;
    case THREAD_PER_SITE_LIKELIHOODS:      
      {
	int 
	  i; 

	/* compute per-site log likelihoods for the sites/partitions 
	   that are handled by this thread */
	perSiteLogLikelihoodsPthreads(localTree, localPr, localTree->lhs, n, tid);

	/* do a parallel gather operation, the threads will write their results 
	   into the global buffer tr->lhs that will then contain all per-site log likelihoods
	   in the proper order 
	*/

	collectDouble(tr->lhs,                localTree->lhs,                  localTree, localPr, n, tid);

      }
      break;
      /* check for errors */
    case THREAD_NEWVIEW_ANCESTRAL:       
      assert(0);
      break; 
    case THREAD_GATHER_ANCESTRAL:
      assert(0); 
      break; 
    case THREAD_EXIT_GRACEFULLY: 
      {
#ifdef MEASURE_TIME_PARALLEL
	double
	  mins[NUM_PAR_JOBS], maxs[NUM_PAR_JOBS]; 
	reduceTimesWorkerRegions(tr, mins, maxs); 
	if(MASTER_P)
	  printParallelTimePerRegion(mins, maxs);
#endif
#ifdef _FINE_GRAIN_MPI
	MPI_Finalize(); 
#endif
	return FALSE; 
      }
      break; 
    default:
      printf("Job %d\n", currentJob);
      assert(0);
    }

#ifdef MEASURE_TIME_PARALLEL 
  timeBuffer[currentJob] += (gettime() - timeForParallelRegion);   
#endif
  
  return TRUE; 
}




/**
   @brief a buziness thread for pthread worker. 
   
   @param tData a struct with basic thread info 
 */ 
void *likelihoodThread(void *tData)
{
  threadData *td = (threadData*)tData;
  tree
    *tr = td->tr;
  partitionList *pr = td->pr;

#ifdef _USE_PTHREADS
  tree *localTree = rax_calloc(1,sizeof(tree )); 
  partitionList *localPr = rax_calloc(1,sizeof(partitionList));

  int
    myCycle = 0;

  const int 
    n = td->tr->numberOfThreads,
    tid = td->threadNumber;

#ifndef _PORTABLE_PTHREADS
  pinToCore(tid);
#endif

  printf("\nThis is RAxML Worker Pthread Number: %d\n", tid);

  while(1)
    {

      while (myCycle == threadJob);
      myCycle = threadJob;

      if ((threadJob >> 16) != THREAD_INIT_PARTITION) {
    	  localPr->perGeneBranchLengths = pr->perGeneBranchLengths;
      	  localPr->numberOfPartitions = pr->numberOfPartitions;
      }
      execFunction(tr, localTree, pr, localPr, tid, n);

      barrierBuffer[tid] = 1;     
    }
#else 
  const int
    n = processes, 
    tid = ((threadData*)tData)->threadNumber;

  printf("\nThis is RAxML Worker Process Number: %d\n", tid);

  while(execFunction(tr, tr, pr, pr, tid,n));
#endif

  return (void*)NULL;
}





/**
   @brief Cleanup step once the master barrier succeeded. 

   This is master specific code called once the barrier is
   passed. Stuff such as reduction operations.  If we execute this
   here, we can keep the code mostly free from parallel -specific
   code.
   
   @param jobType type of parallel region
   @param tr tree 
*/
void masterPostBarrier(int jobType, tree *tr, partitionList *pr)
{
  assert(tr->threadID == 0); 
  
  switch(jobType)
    {
    case THREAD_EVALUATE: 
    case THREAD_OPT_RATE: 
    case THREAD_OPT_ALPHA: 
      {
#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS
	int i,j;
	volatile double partitionResult;	

	for(j = 0; j < pr->numberOfPartitions; j++)
	  {
	    for(i = 0, partitionResult = 0.0; i < tr->numberOfThreads; i++) 
	      partitionResult += globalResult[i * pr->numberOfPartitions+ j];

	    pr->partitionData[j]->partitionLH = partitionResult;
	  }
#endif      
	break; 
      } 
    case THREAD_PER_SITE_LIKELIHOODS:
      {
	int i; 
	/* now just compute the sum over per-site log likelihoods for error checking */      
	double accumulatedPerSiteLikelihood = 0.; 
	for(i = 0; i < tr->originalCrunchedLength; i++)
	  accumulatedPerSiteLikelihood += tr->lhs[i];

	printf("RESULT: %f\t%f", tr->likelihood, accumulatedPerSiteLikelihood); 
	assert(ABS(tr->likelihood - accumulatedPerSiteLikelihood) < 0.00001);
      }
      break;
    } 
}


/**
   @brief a generic master barrier that serves as an entry point for parallel parts of the code.

   @param jobType type of parallel region 
   @param tr tree 
 */ 
void masterBarrier(int jobType, tree *tr, partitionList *pr)
{
#ifdef MEASURE_TIME_PARALLEL
  assert(jobType < NUM_PAR_JOBS); 
  timePerRegion[NUM_PAR_JOBS]  += gettime()- masterTimePerPhase ; 
  masterTimePerPhase = gettime();
#endif

#ifdef _USE_PTHREADS
  const int 
    n = tr->numberOfThreads;

  tr->td[0].functionType = jobType;

  jobCycle = !jobCycle;
  threadJob = (jobType << 16) + jobCycle;

  execFunction(tr, tr, pr, pr, 0, n);

  int 
    i, 
    sum;

  do
    {
      for(i = 1, sum = 1; i < n; i++)
	sum += barrierBuffer[i];
    }
  while(sum < n);  

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
#else 
  tr->td[0].functionType = jobType; 
  execFunction(tr,tr,pr,pr,0,processes);
#endif

  /* code executed by the master, once the barrier is crossed */
  masterPostBarrier(jobType, tr, pr);

#ifdef MEASURE_TIME_PARALLEL
  timePerRegion[jobType] += gettime() - masterTimePerPhase; 
  masterTimePerPhase = gettime();
#endif
}


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

/**
   @brief Tree  initialization function  for workers.

   @param tr tree 
   @param localTree local tree 
   @param tid worker id 
 */ 
/* encapsulated this, s.t. it becomes more clear, that the pthread-master must not execute this */
static void assignAndInitPart1(tree *localTree, tree *tr, partitionList *localPr, partitionList *pr, int *tid)
{
  size_t
    model; 
  int
    totalLength = 0; 

#ifdef _USE_PTHREADS
  localTree->threadID = *tid; 
  printf("my id is %d\n", *tid); 
  assert(localTree != tr);
  localTree->numberOfThreads = tr->numberOfThreads;
#else  /* => MPI */
  *tid = processID; 
  localTree->threadID = processID; 
  tr->numberOfThreads = processes;
#endif

  int bufSize = (8 + pr->numberOfPartitions* 8) * sizeof(int);
  char buf[bufSize], 
    *bufPtr = buf;  
  RECV_BUF(buf, bufSize, MPI_BYTE); 


  ASSIGN_BUF( localTree->useRecom,                  tr->useRecom, int);
  ASSIGN_BUF( localTree->rateHetModel,              tr->rateHetModel, int);
  ASSIGN_BUF( localTree->useMedian,                 tr->useMedian, int); 
  ASSIGN_BUF( localTree->saveMemory,                tr->saveMemory, int);
  ASSIGN_BUF( localTree->maxCategories,             tr->maxCategories, int);
  ASSIGN_BUF( localTree->originalCrunchedLength,    tr->originalCrunchedLength, int);
  ASSIGN_BUF( localTree->mxtips,                    tr->mxtips, int);
  //DIEGO: CHECK THIS
  printf("[DEBUG] ASSIGN PR\n");
  ASSIGN_BUF( localPr->numberOfPartitions,          pr->numberOfPartitions, int);
  ASSIGN_BUF( localPr->perGeneBranchLengths,        pr->perGeneBranchLengths, boolean);
  printf("[DEBUG] ASSIGNED PR\n");

  localTree->td[0].count = 0; 

  if(NOT MASTER_P)
    {
      localTree->lhs                     = (double*)rax_malloc(sizeof(double)   * (size_t)localTree->originalCrunchedLength);     
      localPr->partitionData           = (pInfo**)rax_malloc(NUM_BRANCHES*sizeof(pInfo*));
      for(model = 0; model < (size_t)localPr->numberOfPartitions; model++) {
    	localPr->partitionData[model] = (pInfo*)rax_malloc(sizeof(pInfo));
      }
      localTree->td[0].ti              = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * (size_t)localTree->mxtips);
      localTree->td[0].executeModel    = (boolean *)rax_malloc(sizeof(boolean) * NUM_BRANCHES);
      localTree->td[0].parameterValues = (double *)rax_malloc(sizeof(double) * NUM_BRANCHES);
      localTree->patrat       = (double*)rax_malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);
      localTree->patratStored = (double*)rax_malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);            
    }
  
  for(model = 0; model < (size_t)localPr->numberOfPartitions; model++)
    {
      ASSIGN_BUF(      localPr->partitionData[model]->numberOfCategories,     pr->partitionData[model]->numberOfCategories, int);
      ASSIGN_BUF(      localPr->partitionData[model]->states,                 pr->partitionData[model]->states, int);
      ASSIGN_BUF(      localPr->partitionData[model]->maxTipStates ,          pr->partitionData[model]->maxTipStates, int);
      ASSIGN_BUF(      localPr->partitionData[model]->dataType ,              pr->partitionData[model]->dataType, int);
      ASSIGN_BUF(      localPr->partitionData[model]->protModels ,            pr->partitionData[model]->protModels, int);
      ASSIGN_BUF(      localPr->partitionData[model]->protFreqs ,             pr->partitionData[model]->protFreqs, int);
      ASSIGN_BUF(      localPr->partitionData[model]->lower ,                 pr->partitionData[model]->lower, int);
      ASSIGN_BUF(      localPr->partitionData[model]->upper ,                 pr->partitionData[model]->upper, int);

      localPr->partitionData[model]->partitionLH = 0.0;
      totalLength += (localPr->partitionData[model]->upper -  localPr->partitionData[model]->lower);
    }

  SEND_BUF(buf, bufSize, MPI_BYTE); 

  assert(totalLength == localTree->originalCrunchedLength);

  ASSIGN_DBL(localTree->vectorRecomFraction, tr->vectorRecomFraction); 
}
#endif


/**
   @brief Distribute y-vectors during initialization. 

   @param tr tree 
   @param localTree local tree 
 */ 
void distributeYVectors(tree *localTree, tree *tr, partitionList *localPr, partitionList *pr)
{
  size_t 
    i,
    n = localTree->numberOfThreads,
    globalCounter = 0,
    localCounter = 0,
    model = 0, 
    j; 
  int tid = localTree->threadID; 
  

  /* distribute the y-vectors */
  for(j = 1 ; j <= (size_t)localTree->mxtips; j++)	
    {
#ifdef _FINE_GRAIN_MPI
      unsigned char yBuf[tr->originalCrunchedLength]; 	  
      if(MASTER_P)
	memcpy(yBuf, tr->yVector[j], tr->originalCrunchedLength * sizeof(unsigned char));
      MPI_Bcast(  yBuf, tr->originalCrunchedLength, MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD); 
#endif	  

      for(model = 0, globalCounter = 0; model < (size_t)localPr->numberOfPartitions; model++)
	{
	  if(tr->manyPartitions)
	    {
	      if(isThisMyPartition(localPr, tid, model))
		{
		  assert(localPr->partitionData[model]->upper - localPr->partitionData[model]->lower == localPr->partitionData[model]->width);
		  for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, localCounter++, globalCounter++)
#ifdef _USE_PTHREADS
		    localPr->partitionData[model]->yVector[j][localCounter] = tr->yVector[j][globalCounter];
#else 
		  localPr->partitionData[model]->yVector[j][localCounter] = yBuf[globalCounter];
#endif


		}
	      else
		globalCounter += (localPr->partitionData[model]->upper - localPr->partitionData[model]->lower);
	    }
	  else 
	    {
	      for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, globalCounter++)
		{
		  if(i % (size_t)n == (size_t)tid)
		    {
#ifdef _USE_PTHREADS
		      localPr->partitionData[model]->yVector[j][localCounter] = tr->yVector[j][globalCounter];
#else 
		      localPr->partitionData[model]->yVector[j][localCounter] = yBuf[globalCounter];
#endif
		      ++localCounter; 
		    }
		}	   
	    }
	}
    }
}



/**
   @brief Distribute the weights in the alignment ot workers. 

   @param tr tree 
   @param localTree local tree 
 */ 
void distributeWeights(tree *localTree, tree *tr, partitionList *localPr, partitionList *pr)
{
  int tid = localTree->threadID; 
  int n = localTree->numberOfThreads; 

  size_t     
    globalCounter = 0,
    i,
    localCounter  = 0,
    model; 



  /* distribute the weights  */
#ifdef _FINE_GRAIN_MPI 		/* need to broadcast a few things first */
  if(NOT MASTER_P)
    tr->aliaswgt = rax_malloc(sizeof(int) * tr->originalCrunchedLength); 
  MPI_Bcast(tr->aliaswgt, tr->originalCrunchedLength, MPI_INT, 0, MPI_COMM_WORLD);      
#endif
  for(model = 0, globalCounter = 0; model < (size_t)localPr->numberOfPartitions; model++)
    { 
      if(tr->manyPartitions)
	{
	  if(isThisMyPartition(localPr, tid, model))
	    {
	      assert(localPr->partitionData[model]->upper - localPr->partitionData[model]->lower == localPr->partitionData[model]->width);
	      for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, localCounter++, globalCounter++)
		localPr->partitionData[model]->wgt[localCounter]          = tr->aliaswgt[globalCounter];
	    }
	  else
	    globalCounter += (localPr->partitionData[model]->upper - localPr->partitionData[model]->lower);
	}
      else 
	{ 
	  for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, globalCounter++)
	    {
	      if(i % (size_t)n == (size_t)tid)
		localPr->partitionData[model]->wgt[localCounter++]       = tr->aliaswgt[globalCounter];
	    }	   
	}
    }
}


/**
   @brief Initialize the partitioning scheme (master function).
   @param tr tree 
   @param localTree local tree 
   @param tid worker id    
   @param n number of workers 
 */ 
void initializePartitionsMaster(tree *tr, tree *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
{ 
  size_t
    model;

  treeIsInitialized = TRUE; 

  ASSIGN_INT(localTree->manyPartitions, tr->manyPartitions);
  ASSIGN_INT(localPr->numberOfPartitions, pr->numberOfPartitions);

#ifdef _USE_PTHREADS
  if(MASTER_P)
    globalResult = rax_calloc((size_t) tr->numberOfThreads * (size_t)pr->numberOfPartitions* 2 ,sizeof(double));
  else 
    assignAndInitPart1(localTree, tr, localPr, pr, &tid);
#else 
  globalResult = rax_calloc((size_t) tr->numberOfThreads * (size_t)pr->numberOfPartitions* 2 ,sizeof(double));
  assignAndInitPart1(localTree, tr, localPr, pr, &tid);
  defineTraversalInfoMPI();
#endif

  for(model = 0; model < (size_t)localPr->numberOfPartitions; model++)
    localPr->partitionData[model]->width        = 0;

  if(tr->manyPartitions)    
    {
      multiprocessorScheduling(localTree, localPr, tid);
      computeFractionMany(localPr, tid);
    }
  else
    computeFraction(localPr, tid, n);

  initializePartitionData(localTree, localPr);

  {
    size_t 
      model,  
      i,      
      countOffset,
      myLength = 0;

    for(model = 0; model < (size_t)localPr->numberOfPartitions; model++)
      myLength += localPr->partitionData[model]->width;

    /* assign local memory for storing sequence data */
    
    localTree->y_ptr = (unsigned char *)rax_malloc(myLength * (size_t)(localTree->mxtips) * sizeof(unsigned char));
    assert(localTree->y_ptr != NULL);

    for(i = 0; i < (size_t)localTree->mxtips; i++)
      {
	for(model = 0, countOffset = 0; model < (size_t)localPr->numberOfPartitions; model++)
	  {	    
	    localPr->partitionData[model]->yVector[i+1]   = &localTree->y_ptr[i * myLength + countOffset];
	    countOffset +=  localPr->partitionData[model]->width;
	  }
	assert(countOffset == myLength);
      }

    /* figure in data */

    distributeWeights(localTree, tr, localPr, pr);

    distributeYVectors(localTree, tr, localPr, pr);

  }

  initMemorySavingAndRecom(localTree);
}



