#include "TreeAln.hpp"

#include <limits> 
#include <cassert>

#include "axml.h" 


#if (defined(__SIM_SSE3) && !defined(__AVX))

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#define INTS_PER_VECTOR 4
#define LONG_INTS_PER_VECTOR 2
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128

#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>
#include <pmmintrin.h>

#define INTS_PER_VECTOR 8
#define LONG_INTS_PER_VECTOR 4
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd

#endif


extern const unsigned int mask32[32]; 


// #define REDO_LATER


/************************************************ pop count stuff ***********************************************/

#ifdef REDO_LATER
unsigned int bitcount_32_bit(unsigned int i)
{
  return ((unsigned int) __builtin_popcount(i));
}

/* bit count for 64 bit integers */

inline unsigned int bitcount_64_bit(unsigned long i)
{
  return ((unsigned int) __builtin_popcountl(i));
} 
#endif

/* bit count for 128 bit SSE3 and 256 bit AVX registers */

#if (defined(__SIM_SSE3) || defined(__AVX))
static inline unsigned int vectorPopcount(INT_TYPE v)
{
  unsigned long
    counts[LONG_INTS_PER_VECTOR] __attribute__ ((aligned (BYTE_ALIGNMENT)));

  int    
    i,
    sum = 0;
  
  VECTOR_STORE((CAST)counts, v);

  for(i = 0; i < LONG_INTS_PER_VECTOR; i++)
    sum += __builtin_popcountl(counts[i]);
 	     
  return ((unsigned int)sum);
}
#endif


static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->xPars || (s = s->next)->xPars)
    {
      p->xPars = s->xPars;
      s->xPars = 0;
    }

  assert(p->next->xPars || p->next->next->xPars || p->xPars);
}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, boolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->xPars)
    getxnodeLocal(p);  
  
  if(full)
    {
      if(q->number > maxTips) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->xPars) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->xPars) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}


#if (defined(__SIM_SSE3) || defined(__AVX))


template<unsigned int NUM>
nat addsomeEval(TreeAln &traln, nat model, nat width, nat qNumber, nat pNumber)
{
  nat sum = 0; 
  static INT_TYPE allOne = SET_ALL_BITS_ONE; 

  parsimonyNumber
    *left[NUM],  
    *right[NUM]; 

  auto partition = traln.getPartition(model); 

  // assert(states <= NUM ); 

  for(nat k = 0; k < NUM; k++)
    {
      left[k]  = &(partition->parsVect[(width * NUM * qNumber) + width * k]);
      right[k] = &(partition->parsVect[(width * NUM * pNumber) + width * k]);
    }  
	   
  for(nat i = 0; i < width; i += INTS_PER_VECTOR)
    {                	       
      INT_TYPE      
	l_A,
	v_N = SET_ALL_BITS_ZERO;     
		 
      for(nat j = 0; j < NUM; j++)
	{
	  l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
	  v_N = VECTOR_BIT_OR(l_A, v_N);
	}
		 
      v_N = VECTOR_AND_NOT(v_N, allOne);
      sum += vectorPopcount(v_N); 
    }
  
  return sum; 
}



template<unsigned int NUM>
nat addsome(TreeAln &traln, nat model, nat width, nat qNumber, nat rNumber, nat pNumber)
{
  static INT_TYPE allOne = SET_ALL_BITS_ONE; 
  nat totalScore = 0; 

  auto partition = traln.getPartition(model); 

  parsimonyNumber
    *left[NUM],
    *right[NUM],
    *here[NUM];

  for(nat k = 0; k < NUM; k++)
    {		    
      left[k]  = &(partition->parsVect[(width * NUM * qNumber) + width * k]);
      right[k] = &(partition->parsVect[(width * NUM * rNumber) + width * k]);
      here[k]  = &(partition->parsVect[(width * NUM * pNumber) + width * k]);
    }

  for(nat i = 0; i < width; i += INTS_PER_VECTOR)
    {	 	  
      size_t j;
		    
      INT_TYPE
	s_r, s_l, 
	v_N = SET_ALL_BITS_ZERO,
	l_A[NUM], 
	v_A[NUM];	    	 
		    
      for(j = 0; j < NUM; j++)
	{
	  s_l = VECTOR_LOAD((CAST)(&left[j][i]));
	  s_r = VECTOR_LOAD((CAST)(&right[j][i]));
	  l_A[j] = VECTOR_BIT_AND(s_l, s_r);
	  v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
	  v_N = VECTOR_BIT_OR(v_N, l_A[j]);
	}
		    
      for(j = 0; j < NUM; j++)		    
	VECTOR_STORE((CAST)(&here[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
      v_N = VECTOR_AND_NOT(v_N, allOne);
		    
      totalScore += vectorPopcount(v_N);
    }

  return totalScore; 
} 








static void newviewParsimonyIterativeFast(TreeAln &traln)
{    
  auto tr = traln.getTr(); 

  int 
    *ti = tr->ti,
    count = ti[0],
    index; 

  nat numPart = traln.getNumberOfPartitions();

  for(index = 4; index < count; index += 4)
    {      
      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(nat model = 0; model < numPart; model++)
	{
	  unsigned int
	    totalScore = 0;

	  auto partition = traln.getPartition(model); 

	  size_t
	    states = partition->states,
	    width = partition->parsimonyLength;

	  switch(states)
	    {
	    case 2:       
	      totalScore += addsome<2>(traln, model, width, qNumber, rNumber, pNumber); 
	      break;
	    case 4: 
	      totalScore += addsome<4>(traln, model, width, qNumber, rNumber, pNumber); 
	      break;
	    case 20: 
	      totalScore += addsome<20>(traln, model, width, qNumber, rNumber, pNumber); 
	      break;
	    default:
	      totalScore += addsome<32>(traln, model, width, qNumber, rNumber, pNumber); 
	    }

	  tr->parsimonyScore[pNumber * numPart + model ]
	    = totalScore + tr->parsimonyScore[rNumber * numPart + model] 
	    + tr->parsimonyScore[qNumber * numPart + model]; 	      
	}
    }
}


static void evaluateParsimonyIterativeFast( TreeAln &traln, unsigned int *partitionParsimony, unsigned int *pLengthAtBranch)
{
  auto tr = traln.getTr(); 
  // auto pr = traln.getPartitionsPtr();
  nat numPart = traln.getNumberOfPartitions();

  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(traln);
  
  for(nat i = 0; i < numPart; ++i)
    partitionParsimony[i] = 
      tr->parsimonyScore[pNumber * numPart + i] 
      + tr->parsimonyScore[qNumber * numPart + i]; 

  for(nat model = 0; model < numPart; model++)
    {
      auto partition = traln.getPartition(model);

      size_t
	states = partition->states,
	width = partition->parsimonyLength ;

      nat sum= 0; 

      switch(states)
	{
	case 2:
	  sum += addsomeEval<2>(traln, model, width, qNumber, pNumber); 
	  pLengthAtBranch[model] = sum; 
	  partitionParsimony[model] += sum; 
	  break; 
	case 4:
	  sum += addsomeEval<4>(traln, model, width, qNumber, pNumber); 
	  pLengthAtBranch[model] = sum; 
	  partitionParsimony[model] += sum; 
	  break;
	case 20:
	  sum += addsomeEval<20>(traln, model, width, qNumber, pNumber); 
	  pLengthAtBranch[model] = sum; 
	  partitionParsimony[model] += sum; 
	  break;
	default:
	  sum += addsomeEval<32>(traln, model, width, qNumber, pNumber); 
	  pLengthAtBranch[model] = sum; 
	  partitionParsimony[model] += sum; 
	}
    }
}


#else
#error "compile with either sse3 or avx"
#endif


void evaluateParsimony(TreeAln &traln, nodeptr p, boolean full, unsigned int *partitionParsimony, unsigned int *pLengthAtBranch)
{
  auto tr = traln.getTr(); 

  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->xPars)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->xPars)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }

  ti[0] = counter;
 
  
  evaluateParsimonyIterativeFast(traln ,partitionParsimony, pLengthAtBranch);
}


void newviewParsimony(TreeAln &traln, nodeptr  p)
{     
  auto tr = traln.getTr();

  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(traln);
  }
}


static void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL != 0
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}

static void insertParsimony (TreeAln &traln,  nodeptr p, nodeptr q)
{
  auto tr = traln.getTr(); 

// tree *tr, partitionList *pr
  nodeptr  r;
  
  r = q->back;

  exa_hookupDefault(tr, p->next, q); 
  exa_hookupDefault(tr, p->next->next, r); 

  newviewParsimony(traln, p);
} 


static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  exa_hookupDefault(tr, p,q);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTree (TreeAln &traln, int ip, int iq, int ir)
{    
  // tree *tr, partitionList *pr, 
  auto tr = traln.getTr(); 
  
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  exa_hookupDefault(tr, p,tr->nodep[iq]);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(traln, s,p); 
}


static void makePermutationFast(int *perm, int n, tree *tr)
{    
  int  
    i, 
    j, 
    k;

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {      
      double d =  randum((long int*)(&tr->randomNumberSeed)); // TODO BAD??? 

      k =  (int)((double)(n + 1 - i) * d);
      
      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}


static boolean isInformative(tree *tr, int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const unsigned int
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);	           
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
	informativeCounter++;    
    } 
	  
  if(informativeCounter <= 1)
    return FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(check[j] > 1)
	    return TRUE;
	} 
    }
     
  return FALSE;	     
}


static void determineUninformativeSites(TreeAln &traln, int *informative)
{
  int 
    number = 0; 

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */

  auto tr = traln.getTr(); 

  nat numPart = traln.getNumberOfPartitions();

  for(nat model = 0; model < numPart; model++)
    {
      auto partition = traln.getPartition(model);
      for(nat i = partition->lower; i < nat(partition->upper); i++)
	{
	   if(isInformative(tr, partition->dataType, i))
	     informative[i] = 1;
	   else
	     {
	       informative[i] = 0;
	       number++;
	     }  
	}      
    }

  /* printf("Uninformative Patterns: %d\n", number); */
}


#ifdef REDO_LATER
static void compressDNA(tree *tr, partitionList *pr, int *informative)
{
  size_t
    totalNodes;
    // i,
    // model;
   
  totalNodes = 2 * (size_t)tr->mxtips;

  for(nat model = 0; model < (size_t) pr->numberOfPartitions; model++)
    {
      size_t
	k,
	states = (size_t)pr->partitionData[model]->states,
	compressedEntries,
	compressedEntriesPadded,
	entries = 0, 
	lower = pr->partitionData[model]->lower,
	upper = pr->partitionData[model]->upper;

      parsimonyNumber 
	**compressedTips = (parsimonyNumber **)rax_malloc(states * sizeof(parsimonyNumber*)),
	*compressedValues = (parsimonyNumber *)rax_malloc(states * sizeof(parsimonyNumber));
      
      for(nat i = lower; i < upper; i++)    
	if(informative[i])
	  entries += (size_t)tr->aliaswgt[i];     

      compressedEntries = entries / PCF;

      if(entries % PCF != 0)
	compressedEntries++;

#if (defined(__SIM_SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
	compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
	compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      /* printf("padded is %d\t%d\n", compressedEntriesPadded, compressedEntries);  */

      nat numByte = compressedEntriesPadded * states * totalNodes; 

      /* printf("allocating pars vector of size %d\n", numByte);  */
      pr->partitionData[model]->parsVect = (parsimonyNumber *)rax_malloc_aligned(numByte * sizeof(parsimonyNumber));
     
      for(nat i = 0; i < numByte; i++)      
	pr->partitionData[model]->parsVect[i] = 0;

      for(nat i = 0; i < (size_t)tr->mxtips; i++)
	{
	  size_t
	    w = 0,
	    compressedIndex = 0,
	    compressedCounter = 0,
	    index = 0;

	  for(k = 0; k < states; k++)
	    {
	      compressedTips[k] = &(pr->partitionData[model]->parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
	      compressedValues[k] = 0;
	    }                
	      
	  for(index = lower; index < (size_t)upper; index++)
	    {
	      if(informative[index])
		{
		  const unsigned int 
		    *bitValue = getBitVector(pr->partitionData[model]->dataType);

		  parsimonyNumber 
		    value = bitValue[tr->yVector[i + 1][index]];	  
	      
		  for(w = 0; w < (size_t)tr->aliaswgt[index]; w++)
		    {	   
		      for(k = 0; k < states; k++)
			{
			  if(value & mask32[k])
			    compressedValues[k] |= mask32[compressedCounter];
			}
		     
		      compressedCounter++;
		  
		      if(compressedCounter == PCF)
			{
			  for(k = 0; k < states; k++)
			    {
			      compressedTips[k][compressedIndex] = compressedValues[k];
			      compressedValues[k] = 0;
			    }			 
			  
			  compressedCounter = 0;
			  compressedIndex++;
			}
		    }
		}
	    }
                           
	  for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	    {	
	      for(;compressedCounter < PCF; compressedCounter++)	      
		for(k = 0; k < states; k++)
		  compressedValues[k] |= mask32[compressedCounter];		  
	  
	      for(k = 0; k < states; k++)
		{
		  compressedTips[k][compressedIndex] = compressedValues[k];
		  compressedValues[k] = 0;
		}	      	      
	      
	      compressedCounter = 0;
	    }	 	
	}               

#ifdef DEBUG_PARS
      printf("padded is %d\n", compressedEntriesPadded);
#endif
      pr->partitionData[model]->parsimonyLength = compressedEntriesPadded;

      rax_free(compressedTips);
      rax_free(compressedValues);

      // nat numBytes = totalNodes * pr->partitionData[model]->states * pr->partitionData[model]->parsimonyLength; 

#ifdef DEBUG_PARS
      printf("numBytes=%d\ttotalNodes=%d, tr->partitionData[model].states=%d, parsimonyLength=%d\n",
      	     numBytes,
      	     totalNodes,
      	     pr->partitionData[model]->states,
      	     pr->partitionData[model]->parsimonyLength);

      for(int i = 0; i < numByte; ++i)
      	printf("%u,", pr->partitionData[model]->parsVect[i]);
      printf("\n");
#endif
    }

  tr->parsimonyScore = (unsigned int*)rax_malloc_aligned(sizeof(unsigned int) * totalNodes * pr->numberOfPartitions);  
          
  for(nat i = 0; i < totalNodes * pr->numberOfPartitions; i++) 
    tr->parsimonyScore[i] = 0;
}
#endif


static void stepwiseAddition(TreeAln &traln , nodeptr p, nodeptr q)
{            
// tree *tr, partitionList *pr
  // auto pr = traln.getPartitionsPtr();

  nat numPart = traln.getNumberOfPartitions();
  auto tr = traln.getTr();

  nodeptr 
    r = q->back;

  unsigned int 
    mp;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
  tr->ti[1] = p->number;
  tr->ti[2] = p->back->number;

  nat* partitionParsimony =  (nat*) exa_calloc(numPart,sizeof(nat)) ; 
  nat* pLengthAtBranch = (nat*) exa_calloc(numPart,sizeof(nat)) ; 
  
  evaluateParsimonyIterativeFast(traln, partitionParsimony, pLengthAtBranch);
  
  exa_free(pLengthAtBranch);
  mp = 0; 
  for(nat i = 0; i < numPart; ++i)
    mp += partitionParsimony[i]; 
  exa_free(partitionParsimony); 

  if(mp < tr->bestParsimony)
    {    
      tr->bestParsimony = mp;
      tr->insertNode = q;     
    }
 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {	      
      stepwiseAddition(traln, p, q->next->back);
      stepwiseAddition(traln, p, q->next->next->back);
    }
}


#ifdef REDO_LATER
void allocateParsimonyDataStructures(tree *tr, partitionList *pr)
{
  int 
    i,
    *informative = (int *)rax_malloc(sizeof(int) * (size_t)tr->originalCrunchedLength);
 
  determineUninformativeSites(tr, pr, informative);

  compressDNA(tr, pr, informative);

  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 1; i++)
    {
      nodeptr 
	p = tr->nodep[i];

      p->xPars = 1;		/*  */
      p->next->xPars = 0;
      p->next->next->xPars = 0;
    }

  tr->ti = (int*)rax_malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  rax_free(informative); 
}


void freeParsimonyDataStructures(tree *tr, partitionList *pr)
{
  size_t 
    model;

  rax_free(tr->parsimonyScore);
  
  for(model = 0; model < (size_t) pr->numberOfPartitions; ++model)
    rax_free(pr->partitionData[model]->parsVect);
  
  rax_free(tr->ti);
}
#endif


void makeParsimonyTreeFast(TreeAln &traln )
{   
  auto tr = traln.getTr();

  nodeptr  
    p, 
    f;    

  int 
    nextsp,
    *perm        = (int *)exa_malloc((size_t)(tr->mxtips + 1) * sizeof(int));  

  assert(!tr->constrained);

  makePermutationFast(perm, tr->mxtips, tr);
  
  tr->ntips = 0;    
  
  tr->nextnode = tr->mxtips + 1;       
  
  buildSimpleTree(traln, perm[1], perm[2], perm[3]);
  
  f = tr->start;       
  
  while(tr->ntips < tr->mxtips) 
    {	
      nodeptr q;
      
      tr->bestParsimony = std::numeric_limits<int>::max();
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
        
      if(tr->grouped)
	{
	  int 
	    number = p->back->number;	  	 

	  tr->constraintVector[number] = -9;
	}
          
      stepwiseAddition(traln, q, f->back);
      
      {
	nodeptr	  
	  r = tr->insertNode->back;
	
	int counter = 4;
	
	exa_hookupDefault(tr, q->next,       tr->insertNode);
	exa_hookupDefault(tr, q->next->next, r);

	computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, FALSE);              
	tr->ti[0] = counter;
	
	newviewParsimonyIterativeFast(traln);
      }
    }    
} 

