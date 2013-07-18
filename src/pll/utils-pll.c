/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

/** @file utils.c
 *  
 *  @brief Miscellaneous general utility and helper functions
 */
#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include "cycle.h"


#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
   special bug fix, enforces denormalized numbers to be flushed to zero,
   without this program is a tiny bit faster though.
#include <emmintrin.h> 
#define MM_DAZ_MASK    0x0040
#define MM_DAZ_ON    0x0040
#define MM_DAZ_OFF    0x0000
*/
#endif

#include "axml.h"
#include "phylip_parser/lexer.h"
#include "phylip_parser/phylip.h"
#include "phylip_parser/xalloc.h"
#include "phylip_parser/msa_sites.h"


#include "globalVariables.h"
#include "mem_alloc.h"


extern const unsigned int mask32[32];


/***************** UTILITY FUNCTIONS **************************/

void storeExecuteMaskInTraversalDescriptor(tree *tr, partitionList *pr)
{
  int model;

  for(model = 0; model < pr->numberOfPartitions; model++)
    tr->td[0].executeModel[model] = pr->partitionData[model]->executeModel;

}

void storeValuesInTraversalDescriptor(tree *tr, partitionList *pr, double *value)
{
  int model;

  for(model = 0; model < pr->numberOfPartitions; model++)
    tr->td[0].parameterValues[model] = value[model];
}

#ifdef EXPERIMENTAL
void read_phylip_msa(tree * tr, const char * filename, int format, int type)
{
    size_t
      i, j,
      model;

  struct phylip_data * pd;
  struct msa_sites * ms;
  double **empiricalFrequencies;

  pd = pl_phylip_parse (filename, format);

  ms = construct_msa_sites (pd, SITES_CREATE | SITES_COMPUTE_WEIGHTS);

  free_phylip_struct (pd);
  pd = transpose (ms);
  free_sites_struct (ms);



  tr->mxtips                 = pd->taxa;
  tr->originalCrunchedLength = pd->seqlen;
  pr->numberOfPartitions         = 1;

  setupTree(tr, TRUE);

  tr->gapyness               = 0.03;   /* number of undetermined chars / alignment size */

  /* TODO: The next two lines were commented in model-sep branch */
  tr->aliaswgt = pl_phylip_deldups (&pd);
  tr->originalCrunchedLength = pd->seqlen;
  pr->perGeneBranchLengths = FALSE;

  pl_phylip_subst (pd, DNA_DATA);          /* TODO: Change to reflect the input type */

  tr->rateCategory           =  (int *)    rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));

  tr->patrat                 =  (double *) rax_malloc ((size_t)tr->originalCrunchedLength * sizeof (double));
  tr->patratStored           =  (double *) rax_malloc ((size_t)tr->originalCrunchedLength * sizeof (double));
  tr->lhs                    =  (double *) rax_malloc ((size_t)tr->originalCrunchedLength * sizeof (double));

  tr->executeModel   = (boolean *)rax_malloc(sizeof(boolean) * (size_t)pr->numberOfPartitions);



        
  for(i = 0; i < (size_t)pr->numberOfPartitions; i++)
    tr->executeModel[i] = TRUE;



  /* data structures for convergence criterion need to be initialized after! setupTree */
  if(tr->searchConvergenceCriterion)
  {
    tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
    tr->h = initHashTable(tr->mxtips * 4);
  }

  /* read tip names */
  for(i = 1; i <= (size_t)tr->mxtips; i++)
  {
    tr->nameList[i] = pd->label[i - 1];
  }

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  /* read partition info (boudaries, data type) */
  empiricalFrequencies = (double **)rax_malloc(sizeof(double *) * (size_t)pr->numberOfPartitions);
  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
  {
    int
      len;

    pInfo
      *p = &(pr->partitionData[model]);

    p->states             =  4;   /* TODO: according to the type */
    p->maxTipStates       = 16;   /* TODO: according to the type */
    p->lower              =  0;
    p->upper              = pd->seqlen;
    p->width              = p->upper - p->lower;
    p->dataType           =   DNA_DATA; /* TODO: dna type */
    p->protModels         =   2;
    p->autoProtModels     =   0;
    p->protFreqs          =   0;
    p->nonGTR             =   FALSE;
    p->numberOfCategories =   0;
    
    /* later on if adding secondary structure data

       int    *symmetryVector;
       int    *frequencyGrouping;
       */

    p->partitionName = strdup ("PartName");

//    empiricalFrequencies[model] = (double *)malloc(sizeof(double) * (size_t)pr->partitionData[model]->states);
//    empiricalfrequencies[model][0] = 0.2036082474;
//    empiricalfrequencies[model][1] = 0.2268041237;
//    empiricalfrequencies[model][2] = 0.2731958763;
//    empiricalfrequencies[model][3] = 0.2963917526;
  }
  /* Read all characters from tips */
//  y = (unsigned char *)malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));

  tr->yVector = (char **) rax_malloc(sizeof(char*) * (tr->mxtips+1));
 for (i=0; i < tr->mxtips; ++i)
        tr->yVector[i+1] = pd->seq[i]; //(unsigned char **)malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));
 
 #ifndef _USE_PTHREADS
 #ifndef _FINE_GRAIN_MPI
  //initializePartitionsSequential(tr); 
  initializePartitions (tr, tr, 0, 0);
 #endif
 #endif
}
#endif

/** @brief Read MSA from a file and setup the tree
 *
 *  Reads the MSA from \a filename and constructs the
 *  the tree \a tr and sets up partition and model data
 *
 *  @todo This will be soon replaced by \a read_phylip_msa
 *
 *  @param tr
 *    Pointer to the tree to be set up
 *
 *  @param filename
 *    Filename containing the MSA
 *
 */
void read_msa(tree *tr, partitionList *pr, const char *filename)
  {
    size_t
      i,
      model;

    unsigned char *y;
  double **empiricalFrequencies;

    FILE
      *byteFile = myfopen(filename, "rb");


    /* read the alignment info */
    myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
    myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
    myBinFread(&(pr->numberOfPartitions),         sizeof(int), 1, byteFile);

    /* initialize topology */

    /* Joint branch length estimate is activated by default */
    /*
    if(adef->perGeneBranchLengths)
      tr->numBranches = pr->numberOfPartitions;
    else
      tr->numBranches = 1;
    */
    pr->perGeneBranchLengths = FALSE;
    setupTree(tr, TRUE, pr);
    
    myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);

    /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
       let's not worry about this right now, because it is indeed RAxML-specific */

    tr->aliaswgt                   = (int *)rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
    myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);

    tr->rateCategory    = (int *)    rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
    tr->patrat          = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
    tr->patratStored    = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
    tr->lhs             = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));

    for(i = 0; i < (size_t)pr->numberOfPartitions; i++)
      pr->partitionData[i]->executeModel = TRUE;



    /* data structures for convergence criterion need to be initialized after! setupTree */
    if(tr->searchConvergenceCriterion)
    {
      tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
      tr->h = initHashTable(tr->mxtips * 4);
    }

    /* read tip names */
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

    /* read partition info (boudaries, data type) */
    empiricalFrequencies = (double **)rax_malloc(sizeof(double *) * (size_t)pr->numberOfPartitions);
    for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
    {
      int
        len;

      pInfo
        *p = pr->partitionData[model];

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

      empiricalFrequencies[model] = (double *)rax_malloc(sizeof(double) * (size_t)pr->partitionData[model]->states);
      myBinFread(empiricalFrequencies[model], sizeof(double), pr->partitionData[model]->states, byteFile);
    }
    /* Read all characters from tips */
    y = (unsigned char *)rax_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));

    tr->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength];

    myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

    /* Initialize the model */
    //printf("Here 1!\n");
    initializePartitionsSequential(tr, pr);
    //printf("Here 2!\n");
    initModel(tr, empiricalFrequencies, pr);


    fclose(byteFile);
  }



void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{  
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}

#if 0
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

  ptr = rax_malloc(size);

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
#endif



/* Marked for deletion 
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

*/


void printBothOpen(const char* format, ... )
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

void printResult(tree *tr, partitionList *pr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
  {    
    case TREE_EVALUATION:
      Tree2String(tr->tree_string, tr, pr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, SUMMARIZE_LH, FALSE, FALSE);

      logFile = myfopen(temporaryFileName, "wb");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

      if(adef->perGeneBranchLengths)
        printTreePerGene(tr, pr, adef, temporaryFileName, "wb");
      break;
    case BIG_RAPID_MODE:
      if(finalPrint)
      {
        switch(tr->rateHetModel)
        {
          case GAMMA:
          case GAMMA_I:

            Tree2String(tr->tree_string, tr, pr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint,
                SUMMARIZE_LH, FALSE, FALSE);

            logFile = myfopen(temporaryFileName, "wb");
            fprintf(logFile, "%s", tr->tree_string);
            fclose(logFile);

            if(adef->perGeneBranchLengths)
              printTreePerGene(tr, pr, adef, temporaryFileName, "wb");
            break;
          case CAT:
            /*Tree2String(tr->tree_string, tr, pr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
              NO_BRANCHES, FALSE, FALSE);*/


            Tree2String(tr->tree_string, tr, pr, tr->start->back, TRUE, TRUE, FALSE, FALSE,
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
        Tree2String(tr->tree_string, tr, pr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint,
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



/* Marked for deletion 
boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}
*/

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}

/*
int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}
*/

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}

/*
char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}
*/

const partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
              states    = p->states,
              tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  /*pLength.leftLength = pLength.rightLength = states * states;
    pLength.eignLength = states;
    pLength.evLength   = states * states;
    pLength.eiLength   = states * states;
    pLength.substRatesLength = (states * states - states) / 2;
    pLength.frequenciesLength = states;
    pLength.tipVectorLength   = tipLength * states;
    pLength.symmetryVectorLength = (states * states - states) / 2;
    pLength.frequencyGroupingLength = states;
    pLength.nonGTR = FALSE;*/
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

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

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

/* Marked for deletion 
static int filexists(char *filename)
{
  FILE 
    *fp = fopen(filename,"rb");

  int res; 

  if(fp)
  {
    res = 1;
    fclose(fp);
  }
  else
    res = 0;

  return res;
}
*/


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
  {
    if(fp)
      return fp;
    else
    {	  
      printf("\n Error: the file %s you want to open for reading does not exist, exiting ...\n\n", path);
      exit(-1);
      return (FILE *)NULL;
    }
  }
  else
  {
    if(fp)
      return fp;
    else
    {	 
      printf("\n Error: the file %s you want to open for writing or appending can not be opened [mode: %s], exiting ...\n\n",
          path, mode);
      exit(-1);
      return (FILE *)NULL;
    }
  }


}




/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


/** @brief Check whether a node is a tip.
  * 
  * @param number
  *  Node number to be checked
  *
  * @param maxTips
  *  Number of tips in the tree
  *
  * @return
  *   \b TRUE if tip, \b FALSE otherwise
  */
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

    /* printf("setting x of %d to %d\n", p->number, p->back->number);  */
    p->x = s->x;
    s->x = 0;
  }

  assert(p->x);
}


/** @brief Connect two nodes and assign branch lengths 
  * 
  * Connect the two nodes \a p and \a q in each partition \e i with a branch of
  * length \a z[i]
  *
  * @param p
  *   Node \a p
  * 
  * @param q
  *   Node \a q
  *
  * @param numBranches
  *   Number of partitions
  */
void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

/* connects node p with q and assigns the branch lengths z for the whole vector*/
void hookupFull (nodeptr p, nodeptr q, double *z)
{
  //int i;

  p->back = q;
  q->back = p;

  memcpy(p->z, z, NUM_BRANCHES*sizeof(double) );
  memcpy(q->z, z, NUM_BRANCHES*sizeof(double) );
  //for(i = 0; i < numBranches; i++)
  //  p->z[i] = q->z[i] = z[i];

}

/* connect node p with q and assign the default branch lengths */
void hookupDefault (nodeptr p, nodeptr q)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < NUM_BRANCHES; i++)
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

/* removed the static keyword for using this function in the examples */
boolean setupTree (tree *tr, boolean doInit, partitionList *partitions)
{
  nodeptr  p0, p, q;
  int
    i,
    j;

  int
    tips,
    inter; 

  if(doInit)
    init_default(tr);

  tr->bigCutoff = FALSE;

  tr->maxCategories = MAX(4, tr->categories);

  tips  = (size_t)tr->mxtips;
  inter = (size_t)(tr->mxtips - 1);

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].executeModel = (boolean *)rax_malloc(sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr->td[0].parameterValues = (double *)rax_malloc(sizeof(double) * (size_t)NUM_BRANCHES);

  tr->fracchange = -1.0;

  tr->constraintVector = (int *)rax_malloc((2 * (size_t)tr->mxtips) * sizeof(int));

  tr->nameList = (char **)rax_malloc(sizeof(char *) * (tips + 1));


  p0 = (nodeptr)rax_malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr->nodeBaseAddress = p0;


  tr->nodep = (nodeptr *) rax_malloc((2* (size_t)tr->mxtips) * sizeof(nodeptr));
  assert(tr->nodep);    

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
  {
    p = p0++;

    p->hash   =  KISS32(); /* hast table stuff */
    p->x      =  0;
    p->xBips  = 0;
    p->number =  i;
    p->next   =  p;
    p->back   = (node *)NULL;
    p->bInf   = (branchInfo *)NULL;            
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
      p->bInf   = (branchInfo *)NULL;
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

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionSmoothed[i] = FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;

  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  for (i = 0; i < partitions->numberOfPartitions; i++) {
	partitions->partitionData[i] = (pInfo*)rax_malloc (sizeof(pInfo));
	partitions->partitionData[i]->partitionContribution = -1.0;
	partitions->partitionData[i]->partitionLH = 0.0;
	partitions->partitionData[i]->fracchange = 1.0;
  }

  return TRUE;
}


boolean modelExists(char *model, tree *tr)
{
  /********** BINARY ********************/

  if(strcmp(model, "PSR") == 0)
  {
    tr->rateHetModel = CAT;
    return TRUE;
  }

  if(strcmp(model, "GAMMA") == 0)
  {
    tr->rateHetModel = GAMMA;
    return TRUE;
  }


  return FALSE;
}



/*********************************** *********************************************************/


void init_default(tree *tr)
{

  /*********** tr inits **************/

  tr->numberOfThreads = 1; 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;

  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->saveMemory = FALSE;

  tr->manyPartitions = FALSE;

  tr->startingTree = randomTree;
  tr->randomNumberSeed = 12345;

  tr->categories             = 25;

  tr->grouped = FALSE;
  tr->constrained = FALSE;

  tr->gapyness               = 0.0; 
  tr->useMedian = FALSE;
  /* recom */
  tr->useRecom = FALSE;
  tr->rvec = (recompVectors*)NULL;
  /* recom */

  /********* tr inits end*************/

}








/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/



/* Delete it at some point */
void printLog(tree *tr)
{
  FILE *logFile;
  double t;


  t = gettime() - masterTime;

  logFile = myfopen(logFileName, "ab");

  fprintf(logFile, "%f %f\n", t, tr->likelihood);

  fclose(logFile);


}


void getDataTypeString(tree *tr, pInfo *partitionInfo, char typeOfData[1024])
{
  switch(partitionInfo->dataType)
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





/************************************************************************************/


nodeptr pickRandomSubtree(tree *tr)
{
  nodeptr p;
  do
  {
    int exitDirection = rand() % 3; 
    p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];
    switch(exitDirection)
    {
      case 0:
        break;
      case 1:
        p = p->next;
        break;
      case 2:
        p = p->next->next;
        break;
      default:
        assert(0);
    }
  }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));
  return p;
}
/* small example program that executes ancestral state computations 
   on the entire subtree rooted at p.

   Note that this is a post-order traversal.
*/

  
void computeAllAncestralVectors(nodeptr p, tree *tr, partitionList *pr)
{
  /* if this is not a tip, for which evidently it does not make sense 
     to compute the ancestral sequence because we have the real one ....
  */

  if(!isTip(p->number, tr->mxtips))
    {
      /* descend recursively to compute the ancestral states in the left and right subtrees */

      computeAllAncestralVectors(p->next->back, tr, pr);
      computeAllAncestralVectors(p->next->next->back, tr, pr);
      
      /* then compute the ancestral state at node p */

      newviewGenericAncestral(tr, pr, p);

      /* and print it to terminal, the two booleans that are set to true here 
	 tell the function to print the marginal probabilities as well as 
	 a discrete inner sequence, that is, ACGT etc., always selecting and printing 
	 the state that has the highest probability */

      printAncestralState(p, TRUE, TRUE, tr, pr);
    }
}



void initializePartitionData(tree *localTree, partitionList * localPartitions)
{
  /* in ancestralVectorWidth we store the total length in bytes (!) of 
     one conditional likelihood array !
     we need to know this length such that in the pthreads version the master thread can actually 
     gather the scattered ancestral probabilities from the threads such that they can be printed to screen!
  */

  size_t 
    maxCategories = (size_t)localTree->maxCategories;

  size_t 
    ancestralVectorWidth = 0,
    model; 

  int tid  = localTree->threadID; 

  if(tid > 0)
      localTree->rateCategory    = (int *)    rax_calloc((size_t)localTree->originalCrunchedLength, sizeof(int));	    

  for(model = 0; model < (size_t)localPartitions->numberOfPartitions; model++)
    {
      size_t 
	width = localPartitions->partitionData[model]->width;

      const partitionLengths 
	*pl = getPartitionLengths(localPartitions->partitionData[model]);

      /* 
	 globalScaler needs to be 2 * localTree->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
      */

      localPartitions->partitionData[model]->globalScaler    = (unsigned int *)rax_calloc(2 *(size_t)localTree->mxtips, sizeof(unsigned int));

      localPartitions->partitionData[model]->left              = (double *)rax_malloc_aligned((size_t)pl->leftLength * (maxCategories + 1) * sizeof(double));
      localPartitions->partitionData[model]->right             = (double *)rax_malloc_aligned((size_t)pl->rightLength * (maxCategories + 1) * sizeof(double));
      localPartitions->partitionData[model]->EIGN              = (double*)rax_malloc((size_t)pl->eignLength * sizeof(double));
      localPartitions->partitionData[model]->EV                = (double*)rax_malloc_aligned((size_t)pl->evLength * sizeof(double));
      localPartitions->partitionData[model]->EI                = (double*)rax_malloc((size_t)pl->eiLength * sizeof(double));

      localPartitions->partitionData[model]->substRates        = (double *)rax_malloc((size_t)pl->substRatesLength * sizeof(double));
      localPartitions->partitionData[model]->frequencies       = (double*)rax_malloc((size_t)pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->empiricalFrequencies       = (double*)rax_malloc((size_t)pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->tipVector         = (double *)rax_malloc_aligned((size_t)pl->tipVectorLength * sizeof(double));
      localPartitions->partitionData[model]->symmetryVector    = (int *)rax_malloc((size_t)pl->symmetryVectorLength  * sizeof(int));
      localPartitions->partitionData[model]->frequencyGrouping = (int *)rax_malloc((size_t)pl->frequencyGroupingLength  * sizeof(int));

      localPartitions->partitionData[model]->perSiteRates      = (double *)rax_malloc(sizeof(double) * maxCategories);


      localPartitions->partitionData[model]->nonGTR = FALSE;

      localPartitions->partitionData[model]->gammaRates = (double*)rax_malloc(sizeof(double) * 4);
      localPartitions->partitionData[model]->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char*) * ((size_t)localTree->mxtips + 1));


      localPartitions->partitionData[model]->xVector = (double **)rax_calloc(sizeof(double*), (size_t)localTree->mxtips);

      localPartitions->partitionData[model]->xSpaceVector = (size_t *)rax_calloc((size_t)localTree->mxtips, sizeof(size_t));

      localPartitions->partitionData[model]->sumBuffer = (double *)rax_malloc_aligned(width *
									   (size_t)(localPartitions->partitionData[model]->states) *
									   discreteRateCategories(localTree->rateHetModel) *
									   sizeof(double));


      /* data structure to store the marginal ancestral probabilities in the sequential version or for each thread */

      localPartitions->partitionData[model]->ancestralBuffer = (double *)rax_malloc_aligned(width *
										 (size_t)(localPartitions->partitionData[model]->states) *
										 sizeof(double));

      /* count and accumulate how many bytes we will need for storing a full ancestral vector. for this we addf over the per-partition space requirements in bytes */
      /* ancestralVectorWidth += ((size_t)(pr->partitionData[model]->upper - pr->partitionData[model]->lower) * (size_t)(localPartitions->partitionData[model]->states) * sizeof(double)); */
      ancestralVectorWidth += ((size_t)(localPartitions->partitionData[model]->upper - localPartitions->partitionData[model]->lower) * (size_t)(localPartitions->partitionData[model]->states) * sizeof(double));
      /* :TODO: do we have to use the original tree for that   */

      localPartitions->partitionData[model]->wgt = (int *)rax_malloc_aligned(width * sizeof(int));

      /* rateCategory must be assigned using rax_calloc() at start up there is only one rate category 0 for all sites */

      localPartitions->partitionData[model]->rateCategory = (int *)rax_calloc(width, sizeof(int));

      if(width > 0 && localTree->saveMemory)
	{
	  localPartitions->partitionData[model]->gapVectorLength = ((int)width / 32) + 1;
	  assert(4 == sizeof(unsigned int));
	  localPartitions->partitionData[model]->gapVector = (unsigned int*)rax_calloc((size_t)localPartitions->partitionData[model]->gapVectorLength * 2 * (size_t)localTree->mxtips, sizeof(unsigned int));
	  localPartitions->partitionData[model]->gapColumn = (double *)rax_malloc_aligned(((size_t)localTree->mxtips) *
									       ((size_t)(localPartitions->partitionData[model]->states)) *
									       discreteRateCategories(localTree->rateHetModel) * sizeof(double));
	}
      else
	{
	  localPartitions->partitionData[model]->gapVectorLength = 0;
	  localPartitions->partitionData[model]->gapVector = (unsigned int*)NULL;
	  localPartitions->partitionData[model]->gapColumn = (double*)NULL;
	}              
    }
}

int virtual_width( int n ) {
    const int global_vw = 2;
    return (n+1) / global_vw * global_vw;
}


void initMemorySavingAndRecom(tree *tr, partitionList *pr)
{
  tree
    *localTree = tr; 
  partitionList
    *localPartitions = pr;
  size_t model; 

  /* initialize gap bit vectors at tips when memory saving option is enabled */

  if(localTree->saveMemory)
    {
      for(model = 0; model < (size_t)localPartitions->numberOfPartitions; model++)
	{
	  int        
	    undetermined = getUndetermined(localPartitions->partitionData[model]->dataType);

	  size_t
	    i,
	    j,
	    width =  localPartitions->partitionData[model]->width;

	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(localTree->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(localPartitions->partitionData[model]->yVector[j][i] == undetermined)
		    localPartitions->partitionData[model]->gapVector[localPartitions->partitionData[model]->gapVectorLength * j + i / 32] |= mask32[i % 32];
	    }     
	}
    }
  /* recom */
  if(localTree->useRecom)
    allocRecompVectorsInfo(localTree);
  else
    localTree->rvec = (recompVectors*)NULL;
  /* E recom */
}

double get_branch_length(tree *tr, nodeptr p, int partition_id)
{
  //assert(partition_id < tr->numBranches);
  assert(partition_id < NUM_BRANCHES);
  assert(partition_id >= 0);
  assert(tr->fracchange != -1.0);
  double z = p->z[partition_id];
  if(z < zmin) z = zmin;
  if(z > zmax) z = zmax;
  return (-log(z) * tr->fracchange);
}
void set_branch_length(tree *tr, nodeptr p, int partition_id, double bl)
{
  //assert(partition_id < tr->numBranches);
  assert(partition_id < NUM_BRANCHES);
  assert(partition_id >= 0);
  assert(tr->fracchange != -1.0);
  double z;
  z = exp((-1 * bl)/tr->fracchange);
  if(z < zmin) z = zmin;
  if(z > zmax) z = zmax;
  p->z[partition_id] = z;
}

void initializePartitionsSequential(tree *tr, partitionList *pr)
{ 
  size_t
    model;

  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
    assert(pr->partitionData[model]->width == pr->partitionData[model]->upper - pr->partitionData[model]->lower);

  initializePartitionData(tr, pr);

  /* figure in tip sequence data per-site pattern weights */ 
  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
  {
    size_t
      j;
    size_t lower = pr->partitionData[model]->lower;
    size_t width = pr->partitionData[model]->upper - lower;

    for(j = 1; j <= (size_t)tr->mxtips; j++)
    {
      pr->partitionData[model]->yVector[j] = &(tr->yVector[j][pr->partitionData[model]->lower]);
    }

    memcpy((void*)(&(pr->partitionData[model]->wgt[0])),         (void*)(&(tr->aliaswgt[lower])),      sizeof(int) * width);
  }  

  initMemorySavingAndRecom(tr, pr);
}


/* interface to outside  */
void initializePartitions(tree *tr, tree *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
{
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  initializePartitionsMaster(tr,localTree,pr,localPr,tid,n);
#else
  initializePartitionsSequential(tr, pr);
#endif
}

