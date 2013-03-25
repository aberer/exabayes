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

#ifdef _USE_ZLIB

#include <zlib.h>

#endif

#include <mpi.h>

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
#include "globalVariables.h"





/***************** UTILITY FUNCTIONS **************************/

void storeExecuteMaskInTraversalDescriptor(tree *tr)
{
   int model;
      
   for(model = 0; model < tr->NumberOfModels; model++)
     tr->td[0].executeModel[model] = tr->executeModel[model];
}

void storeValuesInTraversalDescriptor(tree *tr, double *value)
{
   int model;
      
   for(model = 0; model < tr->NumberOfModels; model++)
     tr->td[0].parameterValues[model] = value[model];
}




static void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{  
#ifdef _USE_ZLIB
  int 
    s;
  
   const int 
     max = INT_MAX;
   
   if((size * nmemb) > (size_t)max)
    {
      size_t 	
	toRead = size * nmemb,
	offset = 0;
     
      unsigned char 
	*localPtr = (unsigned char*)ptr;

      size_t 
	rest;

      for(offset = 0; offset < toRead - (size_t)max; offset += (size_t)max)
	{
	  s = gzread(byteFile, (void *)(&localPtr[offset]), max);

	  assert(s == max);      
	}
            
      
      rest = (toRead - offset);

      assert(rest <= (size_t)max);

      s = gzread(byteFile, (void *)(&localPtr[offset]), (int)rest);

      assert(s == (int)rest);
    }
  else    
    {
      s = gzread(byteFile, ptr, (unsigned int)(size * nmemb));

      assert(s == (int)(size * nmemb));
    }

#else
  size_t
    bytes_read;
  
  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
#endif
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

  ptr = malloc(size);
  
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






static void printBoth(FILE *f, const char* format, ... )
{
  if(processID == 0)
    {
      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
    }
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

static void printBothOpenDifferentFile(char *fileName, const char* format, ... )
{
  if(processID == 0)
    {
      FILE 
	*f = myfopen(fileName, "ab");
      
      va_list 
	args;
      
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
            
      fclose(f);
    }
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

static int filexists(char *filename)
{
  FILE 
    *fp;
  
  int 
    res;
  
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

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
	    printf("The file %s ExaML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
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

  
  tr->bigCutoff = FALSE;
  
  tr->maxCategories = MAX(4, tr->categories);
  
  tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);
  
  for(i = 0; i < tr->NumberOfModels; i++)
    tr->partitionContributions[i] = -1.0;
  
  tr->perPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);
  
  
  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->perPartitionLH[i] = 0.0;	    
  
 
  
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

 
 
  
  tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double));
  

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)calloc(tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  
            
  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
  tr->td[0].executeModel = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
  tr->td[0].parameterValues = (double *)malloc(sizeof(double) * tr->NumberOfModels);
  
  for(i = 0; i < tr->NumberOfModels; i++)
    tr->fracchanges[i] = -1.0;
  tr->fracchange = -1.0;
  
  tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));
  
  tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));
   

  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }
  
  tr->nodeBaseAddress = p0;


  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr))))
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

  tr->partitionData = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);

  return TRUE;
}














static void initAdef(analdef *adef)
{   
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = 10;
  adef->bestTrav               = 10;
  adef->initialSet             = FALSE; 
  adef->mode                   = BIG_RAPID_MODE; 
  adef->likelihoodEpsilon      = 0.1;
 
  adef->permuteTreeoptimize    = FALSE; 
  adef->perGeneBranchLengths   = FALSE;  
 
  adef->useCheckpoint          = FALSE;
   
#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
#endif

}




static int modelExists(char *model, tree *tr)
{  
   if(strcmp(model, "PSR\0") == 0)
    {
      tr->rateHetModel = CAT;
      return 1;
    }

  if(strcmp(model, "GAMMA\0") == 0)
    {
      tr->rateHetModel = GAMMA;
      return 1;
    }

  
  return 0;
}



static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }




/*********************************** *********************************************************/


static void printVersionInfo(void)
{
  if(processID == 0)
    printf("\n\nThis is %s version %s released by Alexandros Stamatakis on %s.\n\n",  programName, programVersion, programDate); 
}

static void printMinusFUsage(void)
{
  printf("\n");
 

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("\n");

  printf("               \"-f e\": compute the likelihood of a bunch of trees passed via -t\n");
  printf("                this option will do a quick and dirty optimization without re-optimizng\n");
  printf("                the model parameters for each tree\n");

  printf("\n");

  printf("               \"-f E\": compute the likelihood of a bunch of trees passed via -t\n");
  printf("                this option will do a thorough optimization that re-optimizes\n");
  printf("                the model parameters for each tree\n");

  printf("\n");

  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");
  
  printf("\n");

  printf("              DEFAULT for \"-f\": new rapid hill climbing\n");

  printf("\n");
}


static void printREADME(void)
{
  if(processID == 0)
    {
      printVersionInfo();
      printf("\n");  
      printf("\nTo report bugs use the RAxML google group\n");
      printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
      printf("as well as all error messages printed to screen.\n\n\n");
      
      printf("examl|examl-AVX\n");
      printf("      -s binarySequenceFileName\n");
      printf("      -n outputFileNames\n");
      printf("      -m rateHeterogeneityModel\n");
      printf("      -t userStartingTree|-R binaryCheckpointFile\n");
      printf("      [-a]\n");
      printf("      [-B numberOfMLtreesToSave]\n"); 
      printf("      [-c numberOfCategories]\n");
      printf("      [-D]\n");
      printf("      [-e likelihoodEpsilon] \n");
      printf("      [-f d|e|E|o]\n");    
      printf("      [-h] \n");
      printf("      [-i initialRearrangementSetting] \n");
      printf("      [-M]\n");
      printf("      [-Q]\n");
      printf("      [-S]\n");
      printf("      [-v]\n"); 
      printf("      [-w outputDirectory] \n"); 
      printf("\n");  
      printf("      -a      use the median for the discrete approximation of the GAMMA model of rate heterogeneity\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");
      printf("\n");
      printf("      -B      specify the number of best ML trees to save and print to file\n");
      printf("\n");
      printf("      -c      Specify number of distinct rate catgories for ExaML when modelOfEvolution\n");
      printf("              is set to GTRPSR\n");
      printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
      printf("              categories to accelerate computations. \n");
      printf("\n");
      printf("              DEFAULT: 25\n");
      printf("\n");
      printf("      -D      ML search convergence criterion. This will break off ML searches if the relative \n");
      printf("              Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR cycles\n");
      printf("              is smaller or equal to 1%s. Usage recommended for very large datasets in terms of taxa.\n", "%");
      printf("              On trees with more than 500 taxa this will yield execution time improvements of approximately 50%s\n",  "%");
      printf("              While yielding only slightly worse trees.\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");    
      printf("\n");
      printf("      -e      set model optimization precision in log likelihood units for final\n");
      printf("              optimization of model parameters\n");
      printf("\n");
      printf("              DEFAULT: 0.1 \n"); 
      printf("\n");
      printf("      -f      select algorithm:\n");
      
      printMinusFUsage();
 
      printf("\n");
      printf("      -h      Display this help message.\n");
      printf("\n");  
      printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
      printf("              changes phase\n");
      printf("\n");
      printf("      -m      Model of rate heterogeneity\n");
      printf("\n"); 
      printf("              select \"-m PSR\" for the per-site rate category model (this used to be called CAT in RAxML)\n");
      printf("              select \"-m GAMMA\" for the gamma model of rate heterogeneity with 4 discrete rates\n");
      printf("\n");
      printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");
      printf("              Branch lengths for individual partitions will be printed to separate files\n");
      printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");
      printf("\n");
      printf("      -n      Specifies the name of the output file.\n"); 
      printf("\n");
      printf("      -Q      Enable alternative data/load distribution algorithm for datasets with many partitions\n");
      printf("              In particular under PSR this can lead to parallel performance improvements of up to factor 10!\n");
      printf("\n");
      printf("      -R      read in a binary checkpoint file called ExaML_binaryCheckpoint.RUN_ID_number\n");
      printf("\n");
      printf("      -s      Specify the name of the BINARY alignment data file generated by the parser component\n");
      printf("\n");
      printf("      -S      turn on memory saving option for gappy multi-gene alignments. For large and gappy datasets specify -S to save memory\n");
      printf("              This will produce slightly different likelihood values, may be a bit slower but can reduce memory consumption\n");
      printf("              from 70GB to 19GB on very large and gappy datasets\n");
      printf("\n");
      printf("      -t      Specify a user starting tree file name in Newick format\n");
      printf("\n");
      printf("      -v      Display version information\n");
      printf("\n");
      printf("      -w      FULL (!) path to the directory into which ExaML shall write its output files\n");
      printf("\n");
      printf("              DEFAULT: current directory\n");  
      printf("\n\n\n\n");
    }
}




static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("Error character %c not allowed in run ID\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("Error: please provide a string for the run id after \"-n\" \n");
      assert(0);
    }

}

static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    bad_opt    =FALSE,
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",          
    *optarg,
    model[1024] = "",       
    modelChar;

  double 
    likelihoodEpsilon;
  
  int     
    optind = 1,        
    c,
    nameSet = 0,
    treeSet = 0,   
    modelSet = 0, 
    byteFileSet = 0;


  /*********** tr inits **************/ 
 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;

  tr->manyPartitions = FALSE;

  tr->categories             = 25;

  tr->grouped = FALSE;
  tr->constrained = FALSE;

  tr->gapyness               = 0.0; 
  tr->saveBestTrees          = 0;

  tr->useMedian = FALSE;

  /********* tr inits end*************/




  while(!bad_opt && ((c = mygetopt(argc,argv,"R:B:e:c:f:i:m:t:w:n:s:vhMSDQa", &optind, &optarg))!=-1))
    {
    switch(c)
      {    
      case 'a':
	tr->useMedian = TRUE;
	break;
      case 'B':
	sscanf(optarg,"%d", &(tr->saveBestTrees));
	if(tr->saveBestTrees < 0)
	  {
	    printf("Number of best trees to save must be greater than 0!\n");
	    errorExit(-1);	 
	  }
	break;       
      case 'Q':
	tr->manyPartitions = TRUE;   	
	break;
      case 's':		 	
	strcpy(byteFileName, optarg);	 	
	byteFileSet = TRUE;
	/*printf("%s \n", byteFileName);*/
	break;      
      case 'S':
	tr->saveMemory = TRUE;
	break;
      case 'D':
	tr->searchConvergenceCriterion = TRUE;	
	break;
      case 'R':
	adef->useCheckpoint = TRUE;
	strcpy(binaryCheckpointInputName, optarg);
	break;          
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;                                 
      case 'e':
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;
	break;    
      
      case 'v':
	printVersionInfo();
	errorExit(0);
      
      case 'h':
	printREADME();
	errorExit(0);     
      case 'c':
	sscanf(optarg, "%d", &tr->categories);
	break;     
      case 'f':
	sscanf(optarg, "%c", &modelChar);
	switch(modelChar)
	  {	 
	  case 'e':
	    adef->mode = TREE_EVALUATION;
	    tr->fastTreeEvaluation = TRUE;
	    break;
	  case 'E':
	    adef->mode = TREE_EVALUATION;
	    tr->fastTreeEvaluation = FALSE;
	    break;
	  case 'd':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;	  
	  case 'o':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = FALSE;
	    break;	    	  	  	     
	  default:
	    {
	      if(processID == 0)
		{
		  printf("Error select one of the following algorithms via -f :\n");
		  printMinusFUsage();
		}
	      errorExit(-1);
	    }
	  }
	break;
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;
      case 'n':
        strcpy(run_id,optarg);
	analyzeRunId(run_id);
	nameSet = 1;
        break;
      case 'w':
        strcpy(resultDir, optarg);
	resultDirSet = TRUE;
        break;
      case 't':
	strcpy(tree_file, optarg);       
	treeSet = 1;       
	break;     
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, tr) == 0)
	  {
	    if(processID == 0)
	      {
		printf("Rate heterogeneity Model %s does not exist\n\n", model);               
		printf("For per site rates (called CAT in previous versions) use: PSR\n");	
		printf("For GAMMA use: GAMMA\n");		
	      }
	    errorExit(-1);
	  }
	else
	  modelSet = 1;
	break;     
      default:
	errorExit(-1);
      }
    }



  

  if(!byteFileSet)
    {
      if(processID == 0)
	printf("\nError, you must specify a binary format data file with the \"-s\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    {
      if(processID == 0)
	printf("\nError, you must specify a model of rate heterogeneity with the \"-m\" option\n");
      errorExit(-1);
    }

  if(!nameSet)
    {
      if(processID == 0)
	printf("\nError: please specify a name for this run with -n\n");
      errorExit(-1);
    }

  if(!treeSet && !adef->useCheckpoint)
    {
      if(processID == 0)
	{
	  printf("\nError: please either specify a starting tree for this run with -t\n");
	  printf("or re-start the run from a checkpoint with -R\n");
	}
      
      errorExit(-1);
    }
  
   {

    const 
      char *separator = "/";

    if(resultDirSet)
      {
	char 
	  dir[1024] = "";
	

	if(resultDir[0] != separator[0])
	  strcat(dir, separator);
	
	strcat(dir, resultDir);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	strcpy(workdir, dir);
      }
    else
      {
	char 
	  dir[1024] = "",
	  *result = getcwd(dir, sizeof(dir));
	
	assert(result != (char*)NULL);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	
	strcpy(workdir, dir);		
      }
   }

  return;
}




void errorExit(int e)
{
  MPI_Finalize();

  exit(e);
}



static void makeFileNames(void)
{
  int 
    infoFileExists = 0;
    
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);  
  strcpy(infoFileName,         workdir);
  strcpy(binaryCheckpointName, workdir);
  strcpy(modelFileName, workdir);
   
  strcat(resultFileName,       "ExaML_result.");
  strcat(logFileName,          "ExaML_log.");  
  strcat(infoFileName,         "ExaML_info.");
  strcat(binaryCheckpointName, "ExaML_binaryCheckpoint.");
  strcat(modelFileName,        "ExaML_modelFile.");
  
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);  
  strcat(infoFileName,         run_id); 
  strcat(binaryCheckpointName, run_id);
  strcat(modelFileName,        run_id);

  infoFileExists = filexists(infoFileName);

  if(infoFileExists)
    {
      if(processID == 0)
	{
	  printf("ExaML output files with the run ID <%s> already exist \n", run_id);
	  printf("in directory %s ...... exiting\n", workdir);
	}

      errorExit(-1);	
    }
}




 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{
  if(processID == 0)
    {
      int i, model;
      FILE *infoFile = myfopen(infoFileName, "ab");
      char modelType[128];

      
      if(tr->useMedian)
	strcpy(modelType, "GAMMA with Median");
      else
	strcpy(modelType, "GAMMA");   
     
      printBoth(infoFile, "\n\nThis is %s version %s released by Alexandros Stamatakis in %s.\n\n",  programName, programVersion, programDate);
     
      
      
     
      printBoth(infoFile, "\nAlignment has %zu distinct alignment patterns\n\n",  tr->originalCrunchedLength);
      
     
      
      printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");
      

      switch(adef->mode)
	{	
	case  BIG_RAPID_MODE:	 
	  printBoth(infoFile, "\nExaML rapid hill-climbing mode\n\n");
	  break;
	case TREE_EVALUATION:
	  printBoth(infoFile, "\nExaML %s tree evaluation mode\n\n", (tr->fastTreeEvaluation)?"fast":"slow");
	  break;
	default:
	  assert(0);
	}

     
	  
      if(adef->perGeneBranchLengths)
	printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
      else
	printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);	
	

      
     


      
      
      printBoth(infoFile, "All free model parameters will be estimated by ExaML\n");
      
     
	
      if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
	printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
      else
	{
	  printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", tr->categories);
	  /*
	    if(adef->mode != CLASSIFY_ML)
	    printBoth(infoFile, "Likelihood of final tree will be evaluated and optimized under %s\n\n", modelType);
	  */
	}
      
      /*
	if(adef->mode != CLASSIFY_ML)
	printBoth(infoFile, "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n",
	modelType, adef->likelihoodEpsilon);
      */
    
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  printBoth(infoFile, "Partition: %d\n", model);
	  printBoth(infoFile, "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
	  printBoth(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);
	  
	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      printBoth(infoFile, "DataType: DNA\n");	     
	      printBoth(infoFile, "Substitution Matrix: GTR\n");
	      break;
	    case AA_DATA:
	      assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
	      printBoth(infoFile, "DataType: AA\n");	      
	      printBoth(infoFile, "Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);
	      printBoth(infoFile, "%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");	     
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
	      printBoth(infoFile, "DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
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

      printBoth(infoFile, "ExaML was called as follows:\n\n");
      for(i = 0; i < argc; i++)
	printBoth(infoFile,"%s ", argv[i]);
      printBoth(infoFile,"\n\n\n");

      fclose(infoFile);
    }
}

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
		    NO_BRANCHES, FALSE, FALSE);*/
		  
		  
		  
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
static void printRatesDNA_BIN(int n, double *r, char **names, char *fileName)
{
  int i, j, c;

  for(i = 0, c = 0; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	{
	  if(i == n - 2 && j == n - 1)
	    printBothOpenDifferentFile(fileName, "rate %s <-> %s: %f\n", names[i], names[j], 1.0);
	  else
	    printBothOpenDifferentFile(fileName, "rate %s <-> %s: %f\n", names[i], names[j], r[c]);
	  c++;
	}
    }
}

static void printRatesRest(int n, double *r, char **names, char *fileName)
{
  int i, j, c;

  for(i = 0, c = 0; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	{
	  printBothOpenDifferentFile(fileName, "rate %s <-> %s: %f\n", names[i], names[j], r[c]);
	  c++;
	}
    }
}
static double branchLength(int model, double *z, tree *tr)
{
  double x;
  
  x = z[model];
  assert(x > 0);
  if (x < zmin) 
    x = zmin;  
  
 
  assert(x <= zmax);
  
  if(tr->numBranches == 1)             
    x = -log(x) * tr->fracchange;       
  else
    x = -log(x) * tr->fracchanges[model];

  return x;

}


static double treeLengthRec(nodeptr p, tree *tr, int model)
{  
  double 
    x = branchLength(model, p->z, tr);

  if(isTip(p->number, tr->mxtips))  
    return x;    
  else
    {
      double acc = 0;
      nodeptr q;                
     
      q = p->next;      

      while(q != p)
	{
	  acc += treeLengthRec(q->back, tr, model);
	  q = q->next;
	}

      return acc + x;
    }
}

static double treeLength(tree *tr, int model)
{  
  return treeLengthRec(tr->start->back, tr, model);
}

static void printFreqs(int n, double *f, char **names, char *fileName)
{
  int k;

  for(k = 0; k < n; k++)
    printBothOpenDifferentFile(fileName, "freq pi(%s): %f\n", names[k], f[k]);
}

static void printModelParams(tree *tr, analdef *adef, int treeIteration)
{
  int
    model;

  double
    *f = (double*)NULL,
    *r = (double*)NULL;

  char 
    fileName[2048],
    buf[64];

  strcpy(fileName, modelFileName);
  strcat(fileName, ".");
  sprintf(buf, "%d", treeIteration);
  strcat(fileName, buf);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      double tl;
      char typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);      

      printBothOpenDifferentFile(fileName, "\n\n");

      printBothOpenDifferentFile(fileName, "Model Parameters of Partition %d, Name: %s, Type of Data: %s\n",
				 model, tr->partitionData[model].partitionName, typeOfData);
      
      if(tr->rateHetModel == GAMMA)
	printBothOpenDifferentFile(fileName, "alpha: %f\n", tr->partitionData[model].alpha);
     

      if(adef->perGeneBranchLengths)
	tl = treeLength(tr, model);
      else
	tl = treeLength(tr, 0);

      printBothOpenDifferentFile(fileName, "Tree-Length: %f\n", tl);

      f = tr->partitionData[model].frequencies;
      r = tr->partitionData[model].substRates;

      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  {
	    char *freqNames[20] = {"A", "R", "N ","D", "C", "Q", "E", "G",
				   "H", "I", "L", "K", "M", "F", "P", "S",
				   "T", "W", "Y", "V"};

	    printRatesRest(20, r, freqNames, fileName);
	    printBothOpenDifferentFile(fileName, "\n");
	    printFreqs(20, f, freqNames, fileName);
	  }
	  break;
	case DNA_DATA:
	  {
	    char *freqNames[4] = {"A", "C", "G", "T"};

	    printRatesDNA_BIN(4, r, freqNames, fileName);
	    printBothOpenDifferentFile(fileName, "\n");
	    printFreqs(4, f, freqNames, fileName);
	  }
	  break;
	case BINARY_DATA:
	  {
	    char *freqNames[2] = {"0", "1"};

	    printRatesDNA_BIN(2, r, freqNames, fileName);
	    printBothOpenDifferentFile(fileName, "\n");
	    printFreqs(2, f, freqNames, fileName);
	  }
	  break;
	default:
	  assert(0);
	}

      printBothOpenDifferentFile(fileName, "\n");
    }

  printBothOpenDifferentFile(fileName, "\n");
}


static void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(processID == 0)
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
	case TREE_EVALUATION:	
	  printBothOpen("\n\nOverall Time for evaluating the likelihhod of %d trees: %f secs\n\n", tr->numberOfTrees, t); 
	  printBothOpen("\n\nThe model parameters of the trees have been written to files called %s.i\n", modelFileName);
	  printBothOpen("where i is the number of the tree\n\n");
	  printBothOpen("Note that, in case of a restart from a checkpoint, some tree model files will have been produced by previous runs!\n\n");
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

  for(model = 0; model < tr->NumberOfModels; model++)
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

  for(model = 0; model < tr->NumberOfModels; model++)
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
    
    tr->partitionAssignment = (int *)malloc(tr->NumberOfModels * sizeof(int));
    
  for(model = 0; model < tr->NumberOfModels; model++)
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
    printBothOpen("\nMulti-processor partition data distribution enabled (-Q option)\n");

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
	    *assignments = (int *)calloc(n, sizeof(int));  
	  
	  partitionType 
	    *pt = (partitionType *)malloc(sizeof(partitionType) * p);
	  
	  
	  for(i = 0, k = 0; i < tr->NumberOfModels; i++)
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
	      assert(pt[i].partitionNumber >= 0 && pt[i].partitionNumber < tr->NumberOfModels);
	      tr->partitionAssignment[pt[i].partitionNumber] = minIndex;
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
	  
	  free(assignments);
	  free(pt);
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
  
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    tr->partitionData[model].width        = 0;

  if(tr->manyPartitions)
    {
      multiprocessorScheduling(tr, processID);  
      computeFractionMany(tr, processID);
    }
  else
    computeFraction(tr, processID, processes);
  	   
  maxCategories = tr->maxCategories;

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {                       
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[model])); 

      width = tr->partitionData[model].width;
	
      /* 
	 globalScaler needs to be 2 * tr->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
       */

      tr->partitionData[model].globalScaler    = (unsigned int *)calloc(2 * tr->mxtips, sizeof(unsigned int));  	         

      tr->partitionData[model].left              = (double *)malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].right             = (double *)malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].EIGN              = (double*)malloc(pl->eignLength * sizeof(double));
      tr->partitionData[model].EV                = (double*)malloc_aligned(pl->evLength * sizeof(double));
      tr->partitionData[model].EI                = (double*)malloc(pl->eiLength * sizeof(double));
      
      tr->partitionData[model].substRates        = (double *)malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[model].frequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].empiricalFrequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].tipVector         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));
      tr->partitionData[model].symmetryVector    = (int *)malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[model].frequencyGrouping = (int *)malloc(pl->frequencyGroupingLength  * sizeof(int));
      
      tr->partitionData[model].perSiteRates      = (double *)malloc(sizeof(double) * tr->maxCategories);
            
      tr->partitionData[model].nonGTR = FALSE;            

      tr->partitionData[model].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr->partitionData[model].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

      
      tr->partitionData[model].xVector = (double **)malloc(sizeof(double*) * tr->mxtips);   
      	
      for(j = 0; j < (size_t)tr->mxtips; j++)	        	  	  	  	 
	  tr->partitionData[model].xVector[j]   = (double*)NULL;   

      tr->partitionData[model].xSpaceVector = (size_t *)calloc(tr->mxtips, sizeof(size_t));  

      tr->partitionData[model].sumBuffer = (double *)malloc_aligned(width *
									   (size_t)(tr->partitionData[model].states) *
									   discreteRateCategories(tr->rateHetModel) *
									   sizeof(double));
	    
      tr->partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */

      tr->partitionData[model].rateCategory = (int *)calloc(width, sizeof(int));

      if(width > 0 && tr->saveMemory)
	{
	  tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	    
	  tr->partitionData[model].gapVector = (unsigned int*)calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));	  	    	  	  
	    
	  tr->partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)tr->mxtips) *								      
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

        
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    myLength += tr->partitionData[model].width;         
   
  /* assign local memory for storing sequence data */

  tr->y_ptr = (unsigned char *)malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);
   
  for(i = 0; i < (size_t)tr->mxtips; i++)
    {
      for(model = 0, countOffset = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  tr->partitionData[model].yVector[i+1]   = &tr->y_ptr[i * myLength + countOffset];
	  countOffset +=  tr->partitionData[model].width;
	}
      assert(countOffset == myLength);
    }

  /* figure in data */

  if(tr->manyPartitions)
    {
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
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
      
      for(model = 0, globalCounter = 0; model < (size_t)tr->NumberOfModels; model++)
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
   
  /* set up the averaged frac changes per partition such that no further reading accesses to aliaswgt are necessary
     and we can free the array for the GAMMA model */

  if(tr->NumberOfModels > 1)
    {        
      size_t
	*modelWeights = (size_t *)calloc(tr->NumberOfModels, sizeof(size_t)),
	wgtsum = 0;  
      
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)      
	 {
	   size_t
	     lower = tr->partitionData[model].lower,
	     upper = tr->partitionData[model].upper,
	     i;
	   
	   for(i = lower; i < upper; i++)
	     {
	       modelWeights[model] += (size_t)tr->aliaswgt[i];
	       wgtsum              += (size_t)tr->aliaswgt[i];
	     }
	 }
       
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)      	
	 tr->partitionContributions[model] = ((double)modelWeights[model]) / ((double)wgtsum); 
       
       free(modelWeights);
    }

  if(tr->rateHetModel == GAMMA)
    free(tr->aliaswgt);


  y = (unsigned char *)malloc(sizeof(unsigned char) * tr->originalCrunchedLength);

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength), byteFile);
	
      if(tr->manyPartitions)
	{
	  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
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

	  for(model = 0, globalCounter = 0; model < (size_t)tr->NumberOfModels; model++)
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

  free(y);
    
  /* initialize gap bit vectors at tips when memory saving option is enabled */
  
  if(tr->saveMemory)
    {
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
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


static void initializeTree(tree *tr, analdef *adef)
{
  size_t 
    i,
    model;
  
#ifdef _USE_ZLIB
  gzFile
    byteFile = gzopen(byteFileName, "rb");
#else
  FILE 
    *byteFile = fopen(byteFileName, "rb");
#endif

  double 
    **empiricalFrequencies;	 
  
  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(size_t), 1, byteFile);
  myBinFread(&(tr->NumberOfModels),         sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);
   
  empiricalFrequencies = (double **)malloc(sizeof(double *) * tr->NumberOfModels);
  
  if(adef->perGeneBranchLengths)
    tr->numBranches = tr->NumberOfModels;
  else
    tr->numBranches = 1;
  
  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed ExaML-specific */
  
 
    
  tr->aliaswgt                   = (int *)malloc(tr->originalCrunchedLength * sizeof(int));
  myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       
  
  if(tr->rateHetModel == CAT)
    {
      tr->rateCategory    = (int *)    calloc(tr->originalCrunchedLength, sizeof(int));	    
      tr->patrat          = (double*)  malloc(tr->originalCrunchedLength * sizeof(double));
      tr->patratStored    = (double*)  malloc(tr->originalCrunchedLength * sizeof(double)); 
      tr->lhs             = (double*)  malloc(tr->originalCrunchedLength * sizeof(double)); 
    }
  
  tr->executeModel   = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
  
  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    tr->executeModel[i] = TRUE;
   
  setupTree(tr); 
  
  if(tr->searchConvergenceCriterion && processID == 0)
    {                     
      tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
      tr->h = initHashTable(tr->mxtips * 4);     
    }
  
  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int 
	len;
      
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr->nameList[i] = (char*)malloc(sizeof(char) * len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      addword(tr->nameList[i], tr->nameHash, i);        
    }  
 
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
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
      p->partitionName = (char*)malloc(sizeof(char) * len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);
      
      empiricalFrequencies[model] = (double *)malloc(sizeof(double) * p->states);
      myBinFread(empiricalFrequencies[model], sizeof(double), p->states, byteFile);	   
    }     
  
  initializePartitions(tr, byteFile);
  
  
#ifdef _USE_ZLIB
  gzclose(byteFile);
#else
  fclose(byteFile);
#endif

  initModel(tr, empiricalFrequencies); 
 
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    free(empiricalFrequencies[model]);

  free(empiricalFrequencies);
}


static int getNumberOfTrees(char *fileName, boolean getOffsets, long *treeOffsets)
{
  FILE 
    *f = myfopen(fileName, "r");

  int 
    trees = 0,
    ch;

  if(getOffsets)
    treeOffsets[trees] = 0;

  while((ch = fgetc(f)) != EOF)
    {
      if(ch == ';')
	{
	  trees++;
	  if(getOffsets)	    
	    treeOffsets[trees] = ftell(f) + 1;	 	      	   
	}
    }

  assert(trees > 0);

  fclose(f);

  return trees;
}

static void optimizeTrees(tree *tr, analdef *adef)
{
  long 
    *treeOffsets;

  int 
    i;   

  tr->numberOfTrees = getNumberOfTrees(tree_file, FALSE, (long *)NULL);
  

  if(processID == 0)
    accumulatedTime = 0.0;

  treeOffsets = (long *)malloc(sizeof(long) * (tr->numberOfTrees + 1));

  tr->likelihoods = (double *)malloc(sizeof(double) * tr->numberOfTrees);

  getNumberOfTrees(tree_file, TRUE, treeOffsets);
  
  if(processID == 0)   
    printBothOpen("\n\nFound %d trees to evaluate\n\n", tr->numberOfTrees);
  
  i = 0;

  if(adef->useCheckpoint)
    {      
      restart(tr, adef);       		   	    
	  
      i = ckp.treeIteration;
	       
      if(tr->fastTreeEvaluation && i > 0)	
	treeEvaluate(tr, 2);	
      else
	modOpt(tr, 0.1, adef, i);
      
      tr->likelihoods[i] = tr->likelihood;

      if(processID == 0)
	printModelParams(tr, adef, i);
      
      i++;
    }
       
  for(; i < tr->numberOfTrees; i++)
    {     
      FILE 
	*treeFile = myfopen(tree_file, "rb");
	    
      if(fseek(treeFile, treeOffsets[i], SEEK_SET) != 0)
	assert(0);

      tr->likelihood = unlikely;
   
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE);
               
      fclose(treeFile);
 
      tr->start = tr->nodep[1];     
	  
      if(i > 0)
	resetBranches(tr);
      
      evaluateGeneric(tr, tr->start, TRUE);	
      	  
      if(tr->fastTreeEvaluation && i > 0)
	{
	  ckp.state = MOD_OPT;	  	 

	  ckp.treeIteration = i;
	  
	  writeCheckpoint(tr);

	  treeEvaluate(tr, 2);
	}
      else
	{
	  treeEvaluate(tr, 1);      
	  modOpt(tr, 0.1, adef, i);
	}
      
      tr->likelihoods[i] = tr->likelihood;

      if(processID == 0)
	printModelParams(tr, adef, i);
    }

  if(processID == 0)
    for(i = 0; i < tr->numberOfTrees; i++)
      printBothOpen("Likelihood tree %d: %f \n", i, tr->likelihoods[i]);    
}

int main (int argc, char *argv[])
{ 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  printf("\nThis is ExaML FINE-GRAIN MPI Process Number: %d\n", processID);   
  MPI_Barrier(MPI_COMM_WORLD);
  
  {
    tree  *tr = (tree*)malloc(sizeof(tree));
  
    analdef *adef = (analdef*)malloc(sizeof(analdef));   

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
	printBothOpen("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             
      }  
                         
    /* 
       this will re-start ExaML exactly where it has left off from a checkpoint file,
       while checkpointing is important and has to be implemented for the library we should not worry about this right now 
    */
  
    switch(adef->mode)
      {
      case TREE_EVALUATION:
	optimizeTrees(tr, adef);	
	break;
      case BIG_RAPID_MODE:
       
	if(adef->useCheckpoint)
	  {      
	    /* read checkpoint file */
	    restart(tr, adef);       	
	    
	    /* continue tree search where we left it off */
	    computeBIGRAPID(tr, adef, TRUE); 
	  }
	else
	  {
	    /* not important, only used to keep track of total accumulated exec time 
	       when checkpointing and restarts were used */
	    
	    if(processID == 0)
	      accumulatedTime = 0.0;
	    
	    /* get the starting tree: here we just parse the tree passed via the command line 
	   and do an initial likelihood computation traversal 
	   which we maybe should skip, TODO */
	    
	    getStartingTree(tr);     
	    
	    /* 
	       here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
	       This should basically be the first call to the library that actually computes something :-)
	    */
	    
	    evaluateGeneric(tr, tr->start, TRUE);	
	    
	    /* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */
	    
	    treeEvaluate(tr, 1); 
	    
	    /* now start the ML search algorithm */
      
	    
	    computeBIGRAPID(tr, adef, TRUE); 			     
	  }         
	break;
      default:
	assert(0);
      }
      
    /* print some more nonsense into the ExaML_info file */
  
    if(processID == 0)
      finalizeInfoFile(tr, adef);
  }
  
  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}


