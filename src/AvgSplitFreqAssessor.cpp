#include "AvgSplitFreqAssessor.hpp"


/** 
    @file AvgSplitFreqAssessor.hpp

    @brief Calculates the asdsf. 

    @notice This file is full of hacks, to get this somehow going. 
 */ 


#include <fstream>
#include <ncl/ncl.h>

#include <regex>

#include <algorithm> 

#ifdef _ASDSF_STANDALONE
#define _INCLUDE_DEFINITIONS
#endif
#include "globals.h"
#ifdef _ASDSF_STANDALONE
#undef _INCLUDE_DEFINITIONS
#endif

#include "main-common.h"
#include "globalVariables.h"


static void initializeTreeOnly(int numTax, tree **tre ); 
static stringHashtable *initStringHashTableHere(hashNumberType n); 



#ifdef _ASDSF_STANDALONE

boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


nodeptr findAnyTip(nodeptr p, int numsp)
{ 
  return  isTip(p->number, numsp) ? p : findAnyTip(p->next->back, numsp);
}



void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}



#if HAVE_PLL == 0 
void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}
#else 


void hookupDefault (nodeptr p, nodeptr q)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < NUM_BRANCHES; i++)
    p->z[i] = q->z[i] = defaultz;

}

#endif




void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL == 1 
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}


static hashNumberType  hashString(char *p, hashNumberType tableSize)
{
  hashNumberType h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}

int lookupWord(char *s, stringHashtable *h)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return p->nodeNumber;	  	
    }

  return -1;
}

boolean isTip(int number, int mxtips)
{
  return number <= mxtips; 
}
void printBothOpen( const char* format, ... )
{

}
#endif



#if defined( _ASDSF_STANDALONE) && (HAVE_PLL ==  0)
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
#endif




#if HAVE_PLL == 1 
extern "C"
{
  int processID;   
}
#endif






unsigned int **initBitVectorAsdsf(int mxtips, unsigned int *vectorLength)
{
  unsigned int 
    **bitVectors = (unsigned int **)exa_malloc(sizeof(unsigned int*) * 2 * (size_t)mxtips);
  
  int 
    i;

  if(mxtips % MASK_LENGTH == 0)
    *vectorLength = mxtips / MASK_LENGTH;
  else
    *vectorLength = 1 + (mxtips / MASK_LENGTH); 
  
  for(i = 1; i <= mxtips; i++)
    {
      bitVectors[i] = (unsigned int *)exa_calloc((size_t)(*vectorLength), sizeof(unsigned int));
      assert(bitVectors[i]);
      bitVectors[i][(i - 1) / MASK_LENGTH] |= mask32[(i - 1) % MASK_LENGTH];
    }
  
  for(i = mxtips + 1; i < 2 * mxtips; i++) 
    {
      bitVectors[i] = (unsigned int *)exa_malloc(sizeof(unsigned int) * (size_t)(*vectorLength));
      assert(bitVectors[i]);
    }

  return bitVectors;
}



/**
   @brief checks, if everything is in order with the files and the trees it contains  
*/ 
AvgSplitFreqAssessor::AvgSplitFreqAssessor(vector<string> fileNames, int _start, int _end)
  : start(_start)
  , end(_end)
{
  fillTaxaInfo(fileNames[0]); 

  for(string fn : fileNames)
    {
      if(not fileIsCorrect(fn))
	{
	  cerr << "problem parsing file " << fn << endl;  
	  assert(0); 
	}
    }

  tr = NULL; 
  initializeTreeOnly(taxa.size(), &tr );
  unsigned int vectorLength = 0; 
  tr->bitVectors = initBitVectorAsdsf(tr->mxtips, &vectorLength );
  tr->h = initHashTable(10 * tr->mxtips);
  fns = fileNames; 
}







void AvgSplitFreqAssessor::extractBips()
{
  for (auto filename : fns)
    {
      FILE *fh = fopen(filename.c_str()); 
      int c = 0; 
      while( (c = getc(fh)) != "("); 
      ungetc(c, fh); 

      myTreeReadLen(fh, this->tr, FALSE); 

      Tree2stringNexus(tr->tree_string, tr,  tr->nodep[1]->back, 0 ); 
      cout << tr->tree_string << endl; 

      exit(0); 
      
  //     ifstream infile(filename) ; 
  //     string line; 

  //     // regex rx("\\s+tree gen\\.[0-9]+ = \\[..\\] (.)");
  //     regex rx("\\s+.");
  //     int treeNum = 0;  
  //     while( getline(infile, line))
  // 	{
  // 	  smatch s; 
  // 	  if(regex_match(line,s, rx))
  // 	    {
  // 	      treeNum++;
  // 	      // cout << "the tree  is " << s[0] << endl; 
  // 	    }
  // 	}

  //     cout << "found " << treeNum << " trees" << endl; 

    }
}





std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; 
  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}



/**
   @brief checks, if the number <-> taxon mapping is correct and
   whether there are as many trees as requested for the asdsf
   
   @notice I guess this is not how you write a good parser ;-) 
 */
bool AvgSplitFreqAssessor::fileIsCorrect(string fileName)
{
  string whitespace = " \t"; 
  ifstream infile(fileName); 
  string line ; 
  bool foundTaxaStart = false,
    foundTreeStart = false; 
  int numTrees = 0; 

  while(getline(infile, line))
    {
      string cleanline = trim(line); 

      if(foundTreeStart)	// check number of trees 
	{
	  if(cleanline[cleanline.size()-1] == ';') 
	    ++numTrees; 
	  else 
	    break; 		// DONE  
	}
      else if(foundTaxaStart) 	// check if all taxa are there and have the appropriate number 
	{

	  if(cleanline[cleanline.size()-1] == ';')
	    foundTreeStart = true ; 

	  string cleanerString =  trim(cleanline,  ",;"); 
	  int pos = cleanerString.find_first_of(whitespace, 0); 
	  string num = cleanerString.substr(0, pos),
	    name = cleanerString.substr(pos+1, cleanerString.size()); 

	  int index; 
	  istringstream(num) >> index ; 
	  index--; 
	  assert(name.compare(taxa[index]) == 0); // assert same taxa names 
	  
	} 
      else if(cleanline.compare("translate") == 0)
	foundTaxaStart = true; 
    }

  assert(foundTaxaStart && foundTreeStart); 
  return (numTrees >=  (end - start)); 
}



double AvgSplitFreqAssessor::computeAsdsf()
{
  return 0; 
}


void AvgSplitFreqAssessor::fillTaxaInfo(string fileName)
{
  string whiteSpace = " \t"; 
  ifstream infile(fileName); 
  string line; 
  bool foundStart = false; 
  while(getline(infile, line))
    {      
      string cleanLine = trim(line); 

      if(foundStart)
	{
	  bool abort = false; 
	  if(cleanLine[cleanLine.size()-1] == ';') // we are done 
	    abort = true; 

	  string cleanerString = trim(cleanLine, ",;"); 

	  int pos = cleanerString.find_first_of(whiteSpace, 0);
	  string num =  cleanerString.substr(0, pos),
	    name = cleanerString.substr(pos+1, cleanerString.size()); 	  

	  taxa.push_back(name); 
	  if(abort)
	    break; 
	}
      else if(  cleanLine.compare("translate") == 0  )
	foundStart = true; 
    }  
  assert(foundStart); 
}





#ifdef _ASDSF_STANDALONE 
int main(int argc, char** argv)
{
  if(argc < 4)
    {
      cout << "Usage: " << argv[0] << " start end [file ...]  " << endl << endl
	   << "where  " << endl
	   << "start\t is first tree to include (potentially skipping a burnin) " << endl
	   << "end\t is the last generation to include (file may be truncated)" << endl
	   << "[file ...]\t are various ExaBayes_topology* files. " << endl; 
      exit(0); 
    }

  int start = atoi(argv[1]); 
  int end = atoi(argv[2]); 
  
  char *aFile = argv[3]; 	// TODO 


  vector<string> tmp;
  tmp.push_back(string(aFile)); 
  

  AvgSplitFreqAssessor asdsf(tmp,start,end); 
  asdsf.extractBips();

  return 0; 
}


#endif






// TODO reduce all of this  



// TODO remove at some point =/  
#if HAVE_PLL == 0
static boolean setupTree (tree *tr)
{
  nodeptr  p0, p, q;
  int
    i,
    j,   
    tips,
    inter; 

  int numPartitions = tr->NumberOfModels ;
  
  tr->bigCutoff = FALSE;
  
  tr->maxCategories = MAX(4, tr->categories);
  
  tr->partitionContributions = (double *)exa_malloc(sizeof(double) * numPartitions);
  
  for(i = 0; i < tr->NumberOfModels; i++)
    tr->partitionContributions[i] = -1.0;
  
  tr->perPartitionLH = (double *)exa_malloc(sizeof(double) * numPartitions);
  
  
  for(i = 0; i < numPartitions; i++)    
    tr->perPartitionLH[i] = 0.0;	    
  
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  tr->fracchanges  = (double *)exa_malloc(tr->NumberOfModels * sizeof(double));
  

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)exa_calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)exa_calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)exa_calloc(tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  
            
  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)exa_malloc(sizeof(traversalInfo) * tr->mxtips);
  tr->td[0].executeModel = (boolean *)exa_malloc(sizeof(boolean) * tr->NumberOfModels);
  tr->td[0].parameterValues = (double *)exa_malloc(sizeof(double) * tr->NumberOfModels);
  
  for(i = 0; i < tr->NumberOfModels; i++)
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
  
  tr->nameHash = initStringHashTableHere(10 * tr->mxtips);

  tr->partitionData = (pInfo*)exa_malloc(sizeof(pInfo) * tr->NumberOfModels);

  return TRUE;
}
#else 


static boolean setupTree (tree *tr)
{
  nodeptr  p0, p, q;
  int
    i,
    j;

  int
    tips,
    inter; 

  tr->bigCutoff = FALSE;

  tr->maxCategories = MAX(4, tr->categories);

  tips  = (size_t)tr->mxtips;
  inter = (size_t)(tr->mxtips - 1);

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)exa_calloc((size_t)tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)exa_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)exa_calloc((size_t)tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)exa_malloc(sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].executeModel = (boolean *)exa_malloc(sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr->td[0].parameterValues = (double *)exa_malloc(sizeof(double) * (size_t)NUM_BRANCHES);

  tr->fracchange = -1.0;

  tr->constraintVector = (int *)exa_malloc((2 * (size_t)tr->mxtips) * sizeof(int));

  tr->nameList = (char **)exa_malloc(sizeof(char *) * (tips + 1));


  p0 = (nodeptr)exa_malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr->nodeBaseAddress = p0;


  tr->nodep = (nodeptr *) exa_malloc((2* (size_t)tr->mxtips) * sizeof(nodeptr));
  assert(tr->nodep);    

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
  {
    p = p0++;


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

  tr->nameHash = initStringHashTableHere(10 * tr->mxtips);

  return TRUE;
}


#endif



/** 
    @brief only initializes a raw tree, no partitions or alignment information. 
    
    important: does NOT need a bytefile 
 */ 
static void initializeTreeOnly(int numTax, tree **tre )
{
  *tre = (tree*)exa_calloc(1,sizeof(tree)); 
  tree *tr = *tre; 
  
  tr->mxtips = numTax; 

  setupTree(tr);

}



static stringHashtable *initStringHashTableHere(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)exa_malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)exa_calloc((size_t)tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}
